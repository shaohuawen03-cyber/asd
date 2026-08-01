#!/usr/bin/env python3
"""
将 gmx do_dssp -ssdump 输出的二级结构字符串按时间窗口统计
helix / turn / bend 倾向, 复现论文图4.

用法:
    python3 dssp_bins.py ss_pep.sc output.dat [window_ns]

DSSP 码说明:
    H=alpha螺旋  G=3-10螺旋  I=pi螺旋  -> helix
    T=转角(turn)                        -> turn
    S=弯曲(bend)   B=孤立桥(bridge)    -> bend (近似论文的 bend)
    E=伸展链(strand)                    -> 未计入上述三类
    C=无规卷曲(coil)
"""
import sys
import numpy as np

def main():
    src, dst = sys.argv[1], sys.argv[2]
    window = float(sys.argv[3]) if len(sys.argv) > 3 else 50.0

    times, seqs = [], []
    with open(src) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                t = float(parts[0])
            except ValueError:
                continue
            ss = "".join(parts[1:])
            times.append(t)
            seqs.append(ss)

    times = np.array(times)
    nres = len(seqs[0])
    t0, t1 = times[0], times[-1]
    nwin = int(np.ceil((t1 - t0) / window))

    print(f"帧数: {len(times)}, 残基数: {nres}, 时间范围: {t0}-{t1} ns, 窗口: {window} ns -> {nwin} 窗")

    with open(dst, "w") as out:
        out.write("# time_ns helix_frac turn_frac bend_frac\n")
        for i in range(nwin):
            lo = t0 + i * window
            hi = lo + window
            mask = (times >= lo) & (times < hi)
            if mask.sum() == 0:
                continue
            sub = [seqs[j] for j in range(len(times)) if mask[j]]
            helix = turn = bend = 0
            total = 0
            for s in sub:
                for c in s:
                    total += 1
                    if c in "HGI":
                        helix += 1
                    elif c == "T":
                        turn += 1
                    elif c in "SB":
                        bend += 1
            mid = lo + window / 2
            out.write(f"{mid:.1f} {helix/total:.4f} {turn/total:.4f} {bend/total:.4f}\n")

    print(f"已写入 {dst}")

if __name__ == "__main__":
    main()
