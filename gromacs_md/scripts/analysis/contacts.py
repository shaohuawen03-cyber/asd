#!/usr/bin/env python3
"""
复现论文 3.3 节: 非天然接触 (native / non-native contacts)

定义 (论文 3.3):
  - 天然接触(native): 在参考结构(第一帧)中能够找到的接触
  - 非天然接触: 除天然接触外, 距离在 7 A 以内的原子对接触
  - 接触判据: 两个残基任一一对原子距离 < 7 A 即计为一次接触

用法:
    python3 contacts.py [-t top.tpr] [-f md.xtc] [-a a1-a2] [-p p1-p2]

输出:
  inter_contacts.csv  : 每个肽残基与 AChE 的平均非天然接触数 (图5A)
  intra_contacts.csv  : 肽内部每对残基的平均接触数 (图5B)
  frequent_contacts.tsv : 出现 >10 次的非天然接触对 (表1)
"""
import argparse
from pathlib import Path
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

def select_ache_and_pep(u, arg_a, arg_p):
    if arg_a != "default" and arg_p != "default":
        return u.select_atoms(f"resid {arg_a}"), u.select_atoms(f"resid {arg_p}")
    prot = u.select_atoms("protein")
    # 首选: 按链识别 (与 PDB 残基编号无关)。
    # 本项目 PDB 链 A = 4-542 编号(或 1-542), 肽链 B 单独编号 1-7,
    # 旧版 "resid 531-537" 会误选 AChE 自身的 531-537 号残基!
    pep = u.select_atoms("segid B or chainID B")
    if pep.n_residues >= 2:
        return prot - pep, pep
    # 后退: 残基编号 (旧逻辑, 仅对无 chainID 的体系)
    nres = prot.n_residues
    if nres == 530:  # 单独 AChE 单体对照组 (apo, 无肽)
        return u.select_atoms("protein"), u.select_atoms("protein and resid 531-537")
    if nres == 537:  # 7 肽对接复合物体系 (残基编号 1-537 布局)
        return u.select_atoms("resid 1-530"), u.select_atoms("resid 531-537")
    elif nres == 579:  # 论文 42 肽 Aβ(1-42) 体系
        return u.select_atoms("resid 1-537"), u.select_atoms("resid 538-579")
    else:  # 默认假定最后 7 个残基为小肽
        return u.select_atoms(f"resid 1-{nres-7}"), u.select_atoms(f"resid {nres-6}-{nres}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-t", default="md.tpr")
    ap.add_argument("-f", default="md.xtc")
    ap.add_argument("-a", default="default", help="AChE 残基范围 (默认根据蛋白质长度自适应识别)")
    ap.add_argument("-p", default="default", help="肽残基范围 (默认根据蛋白质长度自适应识别)")
    ap.add_argument("-cut", type=float, default=7.0, help="接触距离截断 (A)")
    ap.add_argument("-freq", type=int, default=10, help="频率阈值 (>N 次计入表1)")
    args = ap.parse_args()

    tpr_file = args.t if Path(args.t).exists() else ("md_0_1.tpr" if Path("md_0_1.tpr").exists() else "md.tpr")
    traj_file = "md_fit.xtc" if Path("md_fit.xtc").exists() else ("md_0_1.xtc" if Path("md_0_1.xtc").exists() else ("md_noPBC.xtc" if Path("md_noPBC.xtc").exists() else args.f))
    u = mda.Universe(tpr_file, traj_file)
    a, p = select_ache_and_pep(u, args.a, args.p)
    if len(p) == 0:
        print(">> [单体蛋白模式] 结构中不存在肽原子，跳过该项非天然接触统计。")
        return
    print(f"AChE 原子: {len(a)}, 肽原子: {len(p)}, 帧数: {u.trajectory.n_frames}")

    # 参考结构 (第一帧) 用于定义天然接触
    u.trajectory[0]
    a_ref = a

    def residue_pair_contact(sel1_res, sel2_res):
        """返回 (i, j) 集合, 残基间任一对原子距离 < cut."""
        pairs = set()
        box = sel2_res.dimensions
        d = distance_array(sel1_res.positions, sel2_res.positions, box=box)
        # d 是 (natom1, natom2) 矩阵
        ri = sel1_res.resids
        rj = sel2_res.resids
        idx = np.argwhere(d < args.cut)
        for i, j in idx:
            pairs.add((ri[i], rj[j]))
        return pairs

    # 天然接触 (第一帧)
    u.trajectory[0]
    native_inter = residue_pair_contact(a, p)
    native_intra = residue_pair_contact(p, p)

    inter_count = np.zeros(u.trajectory.n_frames)          # 每帧肽-AChE 接触残基对总数
    intra_count = np.zeros(u.trajectory.n_frames)
    pep_resids = np.unique(p.resids)
    inter_by_res = np.zeros((u.trajectory.n_frames, len(pep_resids)))
    intra_by_res = np.zeros((u.trajectory.n_frames, len(pep_resids)))

    # 频繁非天然接触统计 (只记录非天然)
    from collections import defaultdict
    inter_freq = defaultdict(int)
    intra_freq = defaultdict(int)

    ridx = {r: k for k, r in enumerate(pep_resids)}

    for f, ts in enumerate(u.trajectory):
        # 使用动态跟进坐标的 AtomGroup (无需每帧重新构造查询)
        pa = a
        pp = p
        inter = residue_pair_contact(pa, pp)
        intra = residue_pair_contact(pp, pp)

        nonnat_inter = inter - native_inter
        nonnat_intra = intra - native_intra

        inter_count[f] = len(inter)
        intra_count[f] = len(intra)

        for (i, j) in nonnat_inter:
            if i in ridx:
                inter_by_res[f, ridx[i]] += 1
            if j in ridx:
                inter_by_res[f, ridx[j]] += 1
            inter_freq[(i, j)] += 1
        for (i, j) in nonnat_intra:
            if i in ridx:
                intra_by_res[f, ridx[i]] += 1
            if j in ridx:
                intra_by_res[f, ridx[j]] += 1
            intra_freq[(i, j)] += 1
        if (f + 1) % 500 == 0:
            print(f"  处理帧 {f+1}/{u.trajectory.n_frames}")

    # ---- 图5A: 肽各残基与 AChE 的平均接触数 ----
    with open("inter_contacts.csv", "w") as fo:
        fo.write("pep_resid,avg_contacts\n")
        for k, r in enumerate(pep_resids):
            fo.write(f"{r},{inter_by_res[:, k].mean():.2f}\n")

    # ---- 图5B: 肽内部各残基平均接触 ----
    with open("intra_contacts.csv", "w") as fo:
        fo.write("pep_resid,avg_contacts\n")
        for k, r in enumerate(pep_resids):
            fo.write(f"{r},{intra_by_res[:, k].mean():.2f}\n")

    # ---- 表1: 最频繁非天然接触 (>10 次, 测试模式时自动适应为 10% 帧数) ----
    min_freq = min(args.freq, max(1, int(u.trajectory.n_frames * 0.1)))
    with open("frequent_contacts.tsv", "w") as fo:
        fo.write("type\tres_i\tres_j\tcount\n")
        for (i, j), c in sorted(inter_freq.items(), key=lambda x: -x[1]):
            if c >= min_freq:
                fo.write(f"inter\t{i}\t{j}\t{c}\n")
        for (i, j), c in sorted(intra_freq.items(), key=lambda x: -x[1]):
            if c >= min_freq:
                fo.write(f"intra\t{i}\t{j}\t{c}\n")

    print("平均 每帧肽-AChE 接触对总数: %.1f" % inter_count.mean())
    print("平均 每帧肽内部接触对总数:   %.1f" % intra_count.mean())
    print("输出: inter_contacts.csv, intra_contacts.csv, frequent_contacts.tsv")

if __name__ == "__main__":
    main()
