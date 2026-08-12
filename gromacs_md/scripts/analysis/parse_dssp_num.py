#!/usr/bin/env python3
"""
自动读取 GROMACS 2025 gmx dssp -num 输出的 ss_pep_num.xvg 文件，
将各结构类型原子数按总残基数转换为占比分布，并分箱统计各时间段的
helix / turn / bend / coil 倾向，输出标准的 ss_pep_bins.dat 数据表。

用法:
    python3 parse_dssp_num.py [ss_pep_num.xvg] [ss_pep_bins.dat] [时间窗口ns, 默认50] [总残基数, 默认7]
"""
import sys
from pathlib import Path
import numpy as np

def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("ss_pep_num.xvg")
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("ss_pep_bins.dat")
    window_ns = float(sys.argv[3]) if len(sys.argv) > 3 else 50.0
    nres = float(sys.argv[4]) if len(sys.argv) > 4 else 7.0

    if not src.exists():
        print(f"!!! 提示: 未找到 {src}，无法执行分箱统计。")
        return

    # 读取 xvg 数据行 (跳过注释 #/@)
    times, struc_c, coil_c, sheet_c, bridge_c, bend_c, turn_c, ahel_c, fhel_c, thel_c = [], [], [], [], [], [], [], [], [], []
    with open(src, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@", ";")):
                continue
            parts = line.split()
            if len(parts) < 10:
                continue
            try:
                # GROMACS gmx dssp -num 列顺序:
                # 0:Time(ps) 1:Structure 2:Coil 3:B-Sheet 4:B-Bridge 5:Bend 6:Turn 7:A-Helix 8:5-Helix 9:3-Helix
                t_ns = float(parts[0]) * 0.001
                times.append(t_ns)
                struc_c.append(float(parts[1]))
                coil_c.append(float(parts[2]))
                sheet_c.append(float(parts[3]))
                bridge_c.append(float(parts[4]))
                bend_c.append(float(parts[5]))
                turn_c.append(float(parts[6]))
                ahel_c.append(float(parts[7]))
                fhel_c.append(float(parts[8]))
                thel_c.append(float(parts[9]))
            except ValueError:
                continue

    if not times:
        print(f"!!! 提示: {src} 中没有有效数据点。")
        return

    times = np.array(times)
    # 计算总数量作为分母，如果和不等于 nres 则采用该帧各类型残基总和
    totals = np.array(coil_c) + np.array(sheet_c) + np.array(bridge_c) + np.array(bend_c) + np.array(turn_c) + np.array(ahel_c) + np.array(fhel_c) + np.array(thel_c)
    totals[totals == 0] = nres

    # 合并为主要分类
    helix_frac = (np.array(ahel_c) + np.array(fhel_c) + np.array(thel_c)) / totals
    turn_frac = np.array(turn_c) / totals
    bend_frac = np.array(bend_c) / totals
    coil_frac = (np.array(coil_c) + np.array(sheet_c) + np.array(bridge_c)) / totals

    t0, t1 = times[0], times[-1]
    nwin = max(1, int(np.ceil((t1 - t0) / window_ns)))

    print(f">> [DSSP -num 转换] 轨迹帧数: {len(times)}, 时间范围: {t0:.2f}-{t1:.2f} ns, 窗口: {window_ns} ns -> {nwin} 窗")

    with open(dst, "w", encoding="utf-8") as out:
        out.write("time_ns helix_frac turn_frac bend_frac coil_frac\n")
        for i in range(nwin):
            lo = t0 + i * window_ns
            hi = lo + window_ns
            mask = (times >= lo) & (times < hi)
            if mask.sum() == 0:
                continue
            h_m = float(helix_frac[mask].mean())
            t_m = float(turn_frac[mask].mean())
            b_m = float(bend_frac[mask].mean())
            c_m = float(coil_frac[mask].mean())
            mid = lo + window_ns / 2.0
            out.write(f"{mid:.3f} {h_m:.4f} {t_m:.4f} {b_m:.4f} {c_m:.4f}\n")

    print(f">> [SAVED] 肽段二级结构分布占比数据已写入: {dst}")

if __name__ == "__main__":
    main()
