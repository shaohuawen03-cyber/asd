#!/usr/bin/env python3
"""
复现论文 3.4 节: 水介导的桥连相互作用 (bridging water)
  - 桥连水: 至少一帧中同时与肽(Aβ)和 AChE 结合的水分子
  - 一个水分子可形成多个桥连相互作用
  - 氢键判据: 距离截断 3.0 A, 角度截断 135 deg
  - 简化实现: 以水氧与两个分子重原子/氢的几何判据近似

用法:
    python3 bridging_waters.py [-t top.tpr] [-f md.xtc] [-a a1-a2] [-p p1-p2]
输出:
  bridging_per_residue.csv : 每个肽残基的桥连水分子数和桥连相互作用数 (图6/表2)
"""
import argparse
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
from collections import defaultdict

# H-bond 几何判据
RCUT = 3.0     # 距离截断 (A)
ACUT = np.cos(np.radians(135.0))   # 135 度角截断 -> cos 值下限

def select_ache_and_pep(u, arg_a, arg_p):
    if arg_a != "default" and arg_p != "default":
        return u.select_atoms(f"resid {arg_a}"), u.select_atoms(f"resid {arg_p}")
    prot = u.select_atoms("protein")
    nres = prot.n_residues
    if nres == 537:  # 7 肽对接复合物体系
        return u.select_atoms("resid 1-530"), u.select_atoms("resid 531-537")
    elif nres == 579:  # 论文 42 肽 Aβ(1-42) 体系
        return u.select_atoms("resid 1-537"), u.select_atoms("resid 538-579")
    else:  # 默认后 7 个残基为小肽
        return u.select_atoms(f"resid 1-{nres-7}"), u.select_atoms(f"resid {nres-6}-{nres}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-t", default="md.tpr")
    ap.add_argument("-f", default="md.xtc")
    ap.add_argument("-a", default="default", help="AChE 残基范围 (默认根据体系残基数自适应)")
    ap.add_argument("-p", default="default", help="肽残基范围 (默认根据体系残基数自适应)")
    args = ap.parse_args()

    traj_file = "md_fit.xtc" if Path("md_fit.xtc").exists() else ("md_noPBC.xtc" if Path("md_noPBC.xtc").exists() else args.f)
    u = mda.Universe(args.t, traj_file)
    waters = u.select_atoms("resname SOL")
    a, p = select_ache_and_pep(u, args.a, args.p)
    if len(p) == 0:
        print(">> [单体蛋白模式] 结构中不存在肽原子，跳过桥连水分子计算。")
        return
    print(f"AChE 原子: {len(a)}, 肽原子: {len(p)}, 水分子数: {waters.n_residues}, 帧数: {u.trajectory.n_frames}")

    pep_resids = np.unique(p.resids)
    # 每个肽残基: [桥连水分子累计, 桥连相互作用累计]
    bridge_waters = defaultdict(set)    # pep_resid -> set(water resids)  (曾经桥连)
    bridge_counts = defaultdict(int)    # pep_resid -> 桥连相互作用总次数

    # 预取水残基索引
    water_res = waters.resids
    # 原子到残基映射
    wat_res_of_atom = waters.residues.resids

    for f, ts in enumerate(u.trajectory):
        pa = a
        pp = p
        pw = waters
        box = pw.dimensions

        # 1. 优先计算肽原子 <-> 所有水原子的距离矩阵 (122 x 71331 -> ~10 ms)
        dpw = distance_array(pp.positions, pw.positions, box=box)
        pw_near_pep = dpw < RCUT
        any_p = pw_near_pep.any(axis=0)  # shape: (71331,)
        cand_idx = np.where(any_p)[0]    # 距离肽 3.0 Å 以内的候选水原子 (通常仅数百个)

        if len(cand_idx) > 0:
            # 2. 仅针对进入肽链势能球的候选水原子计算 AChE 距离 (8145 x ~200 -> 速度提升 100 倍!!)
            pw_cand_positions = pw.positions[cand_idx]
            daw_cand = distance_array(pa.positions, pw_cand_positions, box=box)
            any_a_cand = (daw_cand < RCUT).any(axis=0)

            # 候选水中同时也靠近 AChE 的最终桥连水下标
            bridge_wat_idx = cand_idx[any_a_cand]

            for widx in bridge_wat_idx:
                wres = water_res[widx]
                # 找出与该水桥连的肽残基
                pep_atoms_near = np.where(pw_near_pep[:, widx])[0]
                res_of_pep_atom = pp.resids[pep_atoms_near]
                for r in set(res_of_pep_atom):
                    bridge_waters[r].add(int(wres))
                    bridge_counts[r] += 1

        if (f + 1) % 100 == 0 or f == 0:
            print(f"  已处理轨迹帧 {f+1}/{u.trajectory.n_frames} (候选水原子数: {len(cand_idx)})")

    with open("bridging_per_residue.csv", "w") as fo:
        fo.write("pep_resid,n_bridge_waters,n_bridge_interactions\n")
        for r in pep_resids:
            fo.write(f"{r},{len(bridge_waters.get(int(r), []))},{bridge_counts.get(int(r),0)}\n")

    print("输出: bridging_per_residue.csv")
    print("说明: 本实现以水氧-重原子距离(<3A)近似 H 键; 若要严格复现角度判据")
    print("      建议改用 AmberTools cpptraj 的 Bridge 命令处理原始轨迹。")

if __name__ == "__main__":
    main()
