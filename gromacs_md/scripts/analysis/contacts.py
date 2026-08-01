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
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-t", default="md.tpr")
    ap.add_argument("-f", default="md.xtc")
    ap.add_argument("-a", default="1-537", help="AChE 残基范围")
    ap.add_argument("-p", default="538-579", help="肽残基范围")
    ap.add_argument("-cut", type=float, default=7.0, help="接触距离截断 (A)")
    ap.add_argument("-freq", type=int, default=10, help="频率阈值 (>N 次计入表1)")
    args = ap.parse_args()

    u = mda.Universe(args.t, args.f)
    a = u.select_atoms(f"resid {args.a}")
    p = u.select_atoms(f"resid {args.p}")
    print(f"AChE 原子: {len(a)}, 肽原子: {len(p)}, 帧数: {u.trajectory.n_frames}")

    # 参考结构 (第一帧) 用于定义天然接触
    u.trajectory[0]
    a_ref = u.select_atoms(f"resid {args.a}")

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
        # 重新选择原子 (每帧坐标更新)
        pa = u.select_atoms(f"resid {args.a}")
        pp = u.select_atoms(f"resid {args.p}")
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

    # ---- 表1: 最频繁非天然接触 (>10 次) ----
    with open("frequent_contacts.tsv", "w") as fo:
        fo.write("type\tres_i\tres_j\tcount\n")
        for (i, j), c in sorted(inter_freq.items(), key=lambda x: -x[1]):
            if c >= args.freq:
                fo.write(f"inter\t{i}\t{j}\t{c}\n")
        for (i, j), c in sorted(intra_freq.items(), key=lambda x: -x[1]):
            if c >= args.freq:
                fo.write(f"intra\t{i}\t{j}\t{c}\n")

    print("平均 每帧肽-AChE 接触对总数: %.1f" % inter_count.mean())
    print("平均 每帧肽内部接触对总数:   %.1f" % intra_count.mean())
    print("输出: inter_contacts.csv, intra_contacts.csv, frequent_contacts.tsv")

if __name__ == "__main__":
    main()
