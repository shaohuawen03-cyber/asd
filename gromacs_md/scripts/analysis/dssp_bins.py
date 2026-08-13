#!/usr/bin/env python3
"""
将二级结构字符串序列 (每帧一行) 按时间窗口统计 helix / turn / bend / coil 占比。

兼容两种输入格式:
  1. gmx do_dssp -ssdump 输出:  "时间(ps) 二级结构字符串"
  2. gmx dssp -o 输出:          每帧一行纯二级结构字符串 (无时间列, 用帧序号代替)

DSSP 码归并 (与 parse_dssp_num.py 保持一致):
    H G I -> helix (α/3_10/π 螺旋)
    T     -> turn  (转角)
    S B   -> bend  (弯曲/孤立桥)
    其余 (E P ~ = C) -> coil

用法:
    python3 dssp_bins.py ss_pep.sc output.dat [window_ns, 默认50] [时间单位, ps|ns]
"""
import math
import sys


def main():
    src = sys.argv[1]
    dst = sys.argv[2]
    window = float(sys.argv[3]) if len(sys.argv) > 3 else 50.0
    unit = (sys.argv[4] if len(sys.argv) > 4 else "ps").lower()

    time_scale = 0.001 if unit == "ps" else 1.0

    times, seqs = [], []
    with open(src, encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@", ";")):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    t = float(parts[0]) * time_scale
                    ss = "".join(parts[1:])
                except ValueError:
                    # 纯字符串行 (gmx dssp -o 格式): 用帧序号代替时间
                    ss = "".join(parts)
                    t = float(len(times))
            else:
                ss = parts[0]
                t = float(len(times))
            if not ss:
                continue
            times.append(t)
            seqs.append(ss)

    if len(times) == 0:
        print(">> 提示: 二级结构序列数据为空，跳过统计。")
        return

    nres = len(seqs[0])
    t0, t1 = times[0], times[-1]
    span = max(t1 - t0, 1e-9)
    nwin = max(1, int(math.ceil(span / window)))

    print(f"帧数: {len(times)}, 残基数: {nres}, 时间范围: {t0}-{t1}, "
          f"窗口: {window} -> {nwin} 窗")

    with open(dst, "w") as out:
        out.write("time_ns helix_frac turn_frac bend_frac coil_frac\n")
        for i in range(nwin):
            lo = t0 + i * window
            hi = lo + window
            if i == nwin - 1:
                idx = [j for j, t in enumerate(times) if lo <= t <= t1]
            else:
                idx = [j for j, t in enumerate(times) if lo <= t < hi]
            if not idx:
                continue
            helix = turn = bend = coil = total = 0
            for j in idx:
                for c in seqs[j]:
                    total += 1
                    if c in "HGI":
                        helix += 1
                    elif c == "T":
                        turn += 1
                    elif c in "SB":
                        bend += 1
                    else:
                        coil += 1
            mid = lo + window / 2
            out.write(f"{mid:.4f} {helix/total:.4f} {turn/total:.4f} "
                      f"{bend/total:.4f} {coil/total:.4f}\n")

    print(f"已写入 {dst}")


if __name__ == "__main__":
    main()
