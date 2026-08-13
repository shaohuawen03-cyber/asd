#!/usr/bin/env python3
"""
解析 GROMACS gmx dssp -num 输出的 ss_pep_num.xvg 文件，将各二级结构类型的
残基数按时间窗口统计为 helix / turn / bend / coil 占比，输出 ss_pep_bins.dat。

用法:
    python3 parse_dssp_num.py [ss_pep_num.xvg] [ss_pep_bins.dat] [时间窗口ns, 默认50] [总残基数, 默认7]

★ 重要: GROMACS 2023+ 的 gmx dssp -num 真实输出列为 (枚举顺序固定):
    Time, Loops, Breaks, Bends, Turns, PP-II-Helices, π-Helices,
    3_10-Helices, β-Strands, β-Bridges, α-Helices
    旧版本脚本曾误按 "Time, Structure, Coil, B-Sheet, B-Bridge, Bend, Turn,
    A-Helix, 5-Helix, 3-Helix" 解析, 导致 helix 恒为 0、coil 数值错乱。

本脚本优先读取 xvg 头部 "@ sN legend ..." 图例行, 按名称自适应识别每一列;
对无图例或旧版文件则回退到位置映射, 兼容多种 GROMACS 版本输出。

归并规则 (与 dssp_bins.py 保持一致):
    helix = α-螺旋(H) + 3_10-螺旋(G) + π-螺旋(I)
    turn  = 转角(T)
    bend  = 弯曲(S) + 孤立β桥(B)
    coil  = 无规卷曲(~) + 断裂(=) + β-链(E) + PP-II(P)
"""
import math
import re
import sys
from pathlib import Path

# 10 类二级结构键名 (GROMACS 2023+ gmx dssp 内部枚举顺序)
SS_KEYS = ("loop", "break", "bend", "turn", "ppii", "h5", "h3", "strand", "bridge", "ahelix")

# 无图例时的位置回退映射 (键 = 时间列之外的数据列数)
# 10 数据列 = GROMACS 2023+ gmx dssp -num (Loops, Breaks, Bends, Turns, PP-II, π, 3_10, Strand, Bridge, α)
#  9 数据列 = 旧版假定格式 Structure, Coil, B-Sheet, B-Bridge, Bend, Turn, A-Helix, 5-Helix, 3-Helix
#  8 数据列 = 同上但无 5-Helix
POSITIONAL_MAP = {
    10: ("loop", "break", "bend", "turn", "ppii", "h5", "h3", "strand", "bridge", "ahelix"),
    9: (None, "loop", "strand", "bridge", "bend", "turn", "ahelix", "h5", "h3"),
    8: (None, "loop", "strand", "bridge", "bend", "turn", "ahelix", "h3"),
}


def detect_time_scale(path: Path) -> float:
    """根据 xvg 头部 @ xaxis label 判断时间单位, 返回换算为 ns 的缩放系数"""
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if "@" in line and "xaxis" in line and "label" in line:
                    m = re.search(r'label\s+"([^"]*)"', line)
                    if m:
                        label = m.group(1)
                        if "ns" in label:
                            return 1.0
                        if "ps" in label:
                            return 0.001
                    return 0.001
    except OSError:
        pass
    return 0.001


def classify_legend(name: str):
    """把 xvg 图例字符串映射为二级结构键名, 无法识别时返回 None (该列忽略)"""
    n = re.sub(r"\s+", "", name).lower()
    if not n:
        return None
    if "bridge" in n:
        return "bridge"
    if "strand" in n or "sheet" in n:
        return "strand"
    if "turn" in n:
        return "turn"
    if "bend" in n:
        return "bend"
    if "break" in n:
        return "break"
    if "loop" in n or "coil" in n:
        return "loop"
    if "heli" in n:  # 兼容 Helix / Helices (α/3_10/π/PP-II 复数图例)
        if "pp" in n:
            return "ppii"
        if "3" in n or "10" in n:
            return "h3"
        if "5" in n or "pi" in n or "xp" in n or "π" in n:
            return "h5"
        return "ahelix"
    return None


def read_frames(src: Path):
    """读取 xvg 数据行, 返回 (图例字典, 时间列表(ns), 数据行列表)"""
    legends = {}
    scale = detect_time_scale(src)
    times, rows = [], []
    with open(src, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            if line.startswith("@"):
                m = re.search(r'@\s*s(\d+)\s+legend\s+"([^"]*)"', line)
                if m:
                    legends[int(m.group(1))] = m.group(2)
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                t = float(parts[0]) * scale
                vals = [float(p) for p in parts[1:]]
            except ValueError:
                continue
            if len(vals) < 2:
                continue
            times.append(t)
            rows.append(vals)
    return legends, times, rows


def build_column_map(legends, ncols):
    """优先按图例名称识别每一列, 全部失败时按总列数位置回退"""
    mapping = {}
    if legends:
        for idx in range(1, ncols + 1):
            key = classify_legend(legends.get(idx, ""))
            if key:
                mapping[idx] = key
    if not mapping:
        fallback = POSITIONAL_MAP.get(ncols)
        if fallback:
            mapping = {idx: key for idx, key in enumerate(fallback, start=1) if key}
    return mapping


def main():
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("ss_pep_num.xvg")
    dst = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("ss_pep_bins.dat")
    window_ns = float(sys.argv[3]) if len(sys.argv) > 3 else 50.0
    nres = float(sys.argv[4]) if len(sys.argv) > 4 else 7.0

    if not src.exists():
        print(f"!!! 提示: 未找到 {src}，无法执行分箱统计。")
        return

    legends, times, rows = read_frames(src)
    if not times:
        print(f"!!! 提示: {src} 中没有有效数据点。")
        return

    ncols = len(rows[0])
    colmap = build_column_map(legends, ncols)
    print(f">> [DSSP -num 解析] 共 {ncols} 个数据列, 时间单位缩放 x{detect_time_scale(src)}")
    print(f">> [DSSP -num 解析] 图例: {legends}")
    print(f">> [DSSP -num 解析] 列号->类别映射: {colmap}")
    if not colmap:
        print("!!! 警告: 无法识别任何二级结构列, 输出将为全 0。")

    # 每个时间帧计算四大类占比
    fracs = {k: [] for k in ("helix", "turn", "bend", "coil")}
    for vals in rows:
        counts = {k: 0.0 for k in SS_KEYS}
        for idx, key in colmap.items():
            if idx - 1 < len(vals):
                counts[key] += vals[idx - 1]
        total = sum(counts.values())
        if total <= 0:
            total = nres  # 残基总数兜底 (如旧文件仅有 Structure 计数)
        fracs["helix"].append((counts["ahelix"] + counts["h3"] + counts["h5"]) / total)
        fracs["turn"].append(counts["turn"] / total)
        fracs["bend"].append((counts["bend"] + counts["bridge"]) / total)
        fracs["coil"].append(
            (counts["loop"] + counts["break"] + counts["strand"] + counts["ppii"]) / total
        )

    t0, t1 = times[0], times[-1]
    nwin = max(1, int(math.ceil((t1 - t0) / window_ns)))
    print(f">> 轨迹帧数: {len(times)}, 时间范围: {t0:.2f}-{t1:.2f} ns, "
          f"窗口: {window_ns} ns -> {nwin} 窗")

    with open(dst, "w", encoding="utf-8") as out:
        out.write("time_ns helix_frac turn_frac bend_frac coil_frac\n")
        for i in range(nwin):
            lo = t0 + i * window_ns
            hi = lo + window_ns
            if i == nwin - 1:
                sel = [j for j, t in enumerate(times) if lo <= t <= t1]
            else:
                sel = [j for j, t in enumerate(times) if lo <= t < hi]
            if not sel:
                continue
            mid = lo + window_ns / 2.0
            line = [f"{mid:.4f}"]
            for key in ("helix", "turn", "bend", "coil"):
                mean = sum(fracs[key][j] for j in sel) / len(sel)
                line.append(f"{mean:.4f}")
            out.write(" ".join(line) + "\n")

    print(f">> [SAVED] 肽段二级结构分布占比数据已写入: {dst}")


if __name__ == "__main__":
    main()
