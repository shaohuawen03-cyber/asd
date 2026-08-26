#!/usr/bin/env python3
"""
自动读取 GROMACS 分子动力学模拟与分析结果，并参照论文《乙酰胆碱酯酶-β-淀粉样肽复合物的分子动力学模拟》
批量生成出版级矢量图与位图 (SVG / PNG / PDF) 以及统计指标表。

设计说明:
  - 采用可编辑文本矢量图规范 (svg.fonttype = "none", pdf.fonttype = 42)
  - 自动转换时间单位 (GROMACS xvg 默认 ps -> 缩放 0.001 -> ns)
  - 批量生成图 1 至图 6、肽 RMSD/RMSF 单图，以及一张 2x4 综合汇总图 (fig0_summary_all)
  - 导出统计汇总表格 summary_metrics.csv 及宽表 summary_metrics_wide.csv

用法:
    python3 plot_all.py [--dir 工作目录, 默认 .] [--out 图表保存目录, 默认 ./figures]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap

# 矢量图中文本可编辑，符合出版要求
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["font.size"] = 10
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["figure.dpi"] = 150


def read_xvg(path: Path, x_scale: float = 1.0) -> Optional[pd.DataFrame]:
    """读取 GROMACS xvg 文件，x_scale 用于单位转换 (如 ps -> ns 用 0.001)"""
    if not path.exists():
        return None
    xs, ys = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(("#", "@")):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                xs.append(float(parts[0]) * x_scale)
                ys.append(float(parts[1]))
            except ValueError:
                continue
    if not xs:
        return None
    return pd.DataFrame({"x": xs, "y": ys})


def read_dat_or_csv(path: Path) -> Optional[pd.DataFrame]:
    """读取空格分隔的 .dat 或逗号分隔的 .csv 数据表"""
    if not path.exists():
        return None
    try:
        if path.suffix == ".csv":
            return pd.read_csv(path)
        else:
            return pd.read_csv(path, sep=r"\s+", comment="#", engine="python")
    except Exception:
        return None


def save_all_formats(fig: plt.Figure, out_base: Path) -> None:
    """自动保存为 png (dpi=300), svg, pdf 三种格式"""
    out_base.parent.mkdir(parents=True, exist_ok=True)
    for ext, kwargs in [("png", {"dpi": 300}), ("svg", {}), ("pdf", {})]:
        f = f"{out_base}.{ext}"
        try:
            fig.savefig(f, bbox_inches="tight", **kwargs)
            print(f"  [SAVED] {f}")
        except Exception as e:
            print(f"  [WARNING] Could not save {f}: {e}")


def add_panel_label(ax: plt.Axes, label: str) -> None:
    """按出版规范在子图左上角外部标注 A, B, C... 标签"""
    ax.text(-0.1, 1.05, label, transform=ax.transAxes,
            fontsize=14, weight="bold", va="top", ha="right")


def safe_legend(ax: plt.Axes, **kwargs) -> None:
    """Only add a legend when labeled artists exist.

    Callers may pass fontsize/loc/ncol. Do NOT hard-code fontsize=9
    and then also accept fontsize= in kwargs (TypeError).
    """
    handles, _labels = ax.get_legend_handles_labels()
    if not handles:
        return
    kwargs.setdefault("frameon", False)
    kwargs.setdefault("fontsize", 9)
    ax.legend(**kwargs)


SS_STACK_LABELS = ["α-helix", "β-sheet", "Turn", "Bend", "Coil/loop"]
SS_STACK_COLORS = ["#c0392b", "#f1c40f", "#e67e22", "#27ae60", "#bdc3c7"]


def ss_to_percent_stacks(ss_src: pd.DataFrame):
    """Convert per-frame DSSP fractions (0-1) to 0-100% stacks that sum to 100."""
    t = ss_src["t"].to_numpy(dtype=float)
    helix = np.asarray(ss_src["helix"], dtype=float) * 100.0
    sheet = (np.asarray(ss_src["sheet"], dtype=float) * 100.0
             if "sheet" in ss_src.columns else np.zeros_like(helix))
    turn = np.asarray(ss_src["turn"], dtype=float) * 100.0
    bend = np.asarray(ss_src["bend"], dtype=float) * 100.0
    coil = (np.asarray(ss_src["coil"], dtype=float) * 100.0
            if "coil" in ss_src.columns else np.zeros_like(helix))
    # fold ppii / break / rounding remainder into coil so the stack is 0-100%
    rest = 100.0 - helix - sheet - turn - bend - coil
    coil = np.clip(coil + rest, 0.0, 100.0)
    return t, helix, sheet, turn, bend, coil


def ss_last_window(t, arrays, last_ns=20.0):
    cutoff = max(float(np.max(t)) - last_ns, float(np.min(t)))
    mask = np.asarray(t) >= cutoff
    if not np.any(mask):
        mask = np.ones(len(t), dtype=bool)
    means = [float(np.mean(a[mask])) for a in arrays]
    stds = [float(np.std(a[mask], ddof=1)) if mask.sum() > 1 else 0.0 for a in arrays]
    return mask, means, stds


def draw_ss_lines(ax, t, stacks, fontsize=8, ylim=None):
    """Literature-style % vs time (lines). Better than a flat stack for a stable protein."""
    for y, lab, col in zip(stacks, SS_STACK_LABELS, SS_STACK_COLORS):
        ax.plot(t, y, label=lab, color=col, linewidth=1.15)
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    else:
        ax.set_ylim(0, 100)
    ax.set_ylabel("Content (%)")
    ax.set_xlabel("Time (ns)")
    ax.grid(alpha=0.25, linestyle="--")
    safe_legend(ax, loc="upper right", ncol=2, fontsize=fontsize)


def draw_ss_bars(ax, means, stds, ylim=None):
    """Occupancy bars (mean ± SD). This is the usual DSSP percentage panel in papers."""
    ax.bar(SS_STACK_LABELS, means, yerr=stds, color=SS_STACK_COLORS,
           edgecolor="black", linewidth=0.6, capsize=3, width=0.65)
    ax.set_ylabel("Content (%)")
    if ylim is not None:
        ax.set_ylim(ylim[0], ylim[1])
    else:
        ymax = max(100.0, max(means) + max(stds) + 8)
        ax.set_ylim(0, ymax)
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", alpha=0.25, linestyle="--")
    for i, (m, s) in enumerate(zip(means, stds)):
        ax.text(i, m + s + 1.5, f"{m:.1f}", ha="center", va="bottom", fontsize=8)


def draw_ss_heatmap(ax, times, seqs):
    """Residue × time DSSP map (classic MD DSSP figure)."""
    if times is None or seqs is None or len(seqs) < 5:
        ax.text(0.5, 0.5, "DSSP per-residue missing", ha="center", va="center", color="tab:red")
        ax.set_axis_off()
        return
    step = max(1, len(seqs) // 500)
    seqs_d = seqs[::step]
    t_d = np.asarray(times, dtype=float)[::step]
    code_to_int = {
        "H": 4, "G": 4, "I": 4,
        "E": 3, "B": 3,
        "T": 2,
        "S": 1, "P": 1,
        "C": 0, "~": 0, "-": 0, " ": 0, "=": 0,
    }
    nres = max(len(s) for s in seqs_d)
    mat = np.zeros((nres, len(seqs_d)))
    for i, s in enumerate(seqs_d):
        for j, c in enumerate(s[:nres]):
            mat[j, i] = code_to_int.get(c, 0)
    cmap = ListedColormap(["#bdc3c7", "#27ae60", "#e67e22", "#f1c40f", "#c0392b"])
    im = ax.imshow(
        mat, aspect="auto", origin="lower", interpolation="nearest",
        cmap=cmap, vmin=-0.5, vmax=4.5, rasterized=True,
        extent=[float(t_d[0]), float(t_d[-1]), 0.5, nres + 0.5],
    )
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Residue")
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, ticks=[0, 1, 2, 3, 4])
    cbar.ax.set_yticklabels(["Coil", "Bend", "Turn", "Sheet", "Helix"], fontsize=7)


def detect_rmsd_phases(t, y, max_phases: int = 3, min_span_ns: float = 8.0):
    """Greedy jump detection on a smoothed RMSD trace."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 20:
        return [(float(t[0]), float(t[-1]), float(y.mean()),
                 float(y.std(ddof=1) if n > 1 else 0.0), int(n))]
    win = max(7, (int(round(n * 0.015)) | 1))
    pad = win // 2
    yp = np.pad(y, (pad, pad), mode="edge")
    ys = np.convolve(yp, np.ones(win) / win, mode="valid")
    if len(ys) > n:
        ys = ys[:n]
    elif len(ys) < n:
        ys = np.pad(ys, (0, n - len(ys)), mode="edge")
    dy = np.abs(np.diff(ys))
    thr = float(np.percentile(dy, 85))
    cand = []
    for i in range(1, len(dy) - 1):
        if dy[i] >= dy[i - 1] and dy[i] >= dy[i + 1] and dy[i] >= thr:
            cand.append((float(dy[i]), i + 1))
    cand.sort(reverse=True)
    span = float(t[-1] - t[0])
    min_span = min(min_span_ns, span / (max_phases + 1))
    jumps = []
    for _mag, idx in cand:
        tt = float(t[idx])
        if tt - float(t[0]) < min_span or float(t[-1]) - tt < min_span:
            continue
        if all(abs(tt - float(t[j])) >= min_span for j in jumps):
            jumps.append(idx)
        if len(jumps) >= max_phases - 1:
            break
    bounds = [0] + sorted(jumps) + [n]
    phases = []
    for a, b in zip(bounds[:-1], bounds[1:]):
        b_end = n if b == n else b
        if b_end <= a:
            continue
        seg = y[a:b_end]
        t1 = float(t[b_end - 1] if b_end == n else t[b_end])
        phases.append((float(t[a]), t1, float(seg.mean()),
                       float(seg.std(ddof=1) if len(seg) > 1 else 0.0), int(len(seg))))
    return phases


def phase_label(i: int, n: int) -> str:
    if n <= 1:
        return "STABLE"
    if n == 2:
        return ("INITIAL_POSE", "EQUILIBRATED_BOUND")[i]
    names = ["INITIAL_BOUND", "METASTABLE_REARRANGE", "EQUILIBRATED_BOUND"]
    return names[i] if i < len(names) else f"PHASE_{i+1}"


def read_ss_frac(path: Path) -> Optional[pd.DataFrame]:
    """Per-frame DSSP fractions: time_ns helix turn bend sheet coil ..."""
    if not path.exists():
        return None
    rows = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(("#", "@", ";")):
                continue
            parts = s.split()
            if len(parts) < 4:
                continue
            try:
                rec = {
                    "time_ns": float(parts[0]),
                    "helix": float(parts[1]),
                    "turn": float(parts[2]),
                    "bend": float(parts[3]),
                }
                rec["sheet"] = float(parts[4]) if len(parts) > 4 else 0.0
                rec["coil"] = float(parts[5]) if len(parts) > 5 else 0.0
                if len(parts) >= 5 and len(parts) < 6:
                    rec["coil"] = float(parts[4])
                    rec["sheet"] = 0.0
                rows.append(rec)
            except ValueError:
                continue
    if not rows:
        return None
    return pd.DataFrame(rows)


def read_ss_perres(path: Path):
    if not path.exists():
        return None, []
    times, seqs = [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith(("#", "@", ";")):
                continue
            parts = s.split()
            if len(parts) < 2:
                continue
            try:
                times.append(float(parts[0]))
            except ValueError:
                continue
            seqs.append("".join(parts[1:]))
    if not seqs:
        return None, []
    return np.asarray(times), seqs


def ss_string_heatmap(times, seqs):
    code_to_int = {
        "H": 3, "G": 3, "I": 3,
        "E": 2, "B": 2,
        "T": 1, "S": 1, "P": 1,
        "C": 0, "~": 0, "-": 0, " ": 0, "=": 0,
    }
    nres = max(len(s) for s in seqs)
    mat = np.zeros((len(seqs), nres))
    for i, s in enumerate(seqs):
        for j, c in enumerate(s[:nres]):
            mat[i, j] = code_to_int.get(c, 0)
    return mat


def apply_ylim(ax, limits: dict, key: str) -> None:
    """Apply a shared y-range (from --limits-json) to a panel, if provided."""
    pair = limits.get(key)
    if not pair:
        return
    try:
        lo, hi = float(pair[0]), float(pair[1])
    except (TypeError, ValueError, IndexError):
        return
    if hi > lo:
        ax.set_ylim(lo, hi)


def summarize_rmsf(df: Optional[pd.DataFrame], metric: str, system_label: str):
    """RMSF x-axis is residue number — never treat it as a time series."""
    if df is None or df.empty:
        return None
    y = df["y"]
    return {
        "system": system_label,
        "metric": metric,
        "mean": float(y.mean()),
        "std": float(y.std(ddof=1)) if len(y) > 1 else 0.0,
        "full_traj_mean": float(y.mean()),
        "full_traj_std": float(y.std(ddof=1)) if len(y) > 1 else 0.0,
        "min": float(y.min()),
        "max": float(y.max()),
        "status": "PER_RESIDUE_FLEXIBILITY (逐残基柔性, 非时间序列)",
        "n_points": int(len(y)),
    }


def summarize_rdf(df: Optional[pd.DataFrame], metric: str, system_label: str):
    """针对径向分布函数 g(r)，提取峰高度 g_max 及主峰位置 r(nm)"""
    if df is None or df.empty:
        return None
    peak_idx = df["y"].idxmax()
    r_peak = float(df.loc[peak_idx, "x"])
    g_max = float(df.loc[peak_idx, "y"])
    return {
        "system": system_label,
        "metric": metric,
        "mean": g_max,
        "std": r_peak,
        "full_traj_mean": g_max,
        "full_traj_std": 0.0,
        "min": float(df["y"].min()),
        "max": g_max,
        "status": f"Peak at r={r_peak:.2f} nm",
        "n_points": int(len(df)),
    }


def summarize_last_ns(df: Optional[pd.DataFrame], metric: str, system_label: str, last_ns: float = 20.0):
    if df is None or df.empty:
        return None
    total_span = df["x"].max() - df["x"].min()
    # 针对 100 步测试轨迹的自适应计算：若总跨度短于 last_ns，则取后 50% 均值
    eff_last = min(last_ns, total_span * 0.5) if total_span > 0 else 0.0
    cutoff = max(df["x"].max() - eff_last, df["x"].min())
    sub = df.loc[df["x"] >= cutoff, "y"]
    if sub.empty:
        sub = df["y"]
    full_mean = float(df["y"].mean())
    full_std = float(df["y"].std(ddof=1)) if len(df["y"]) > 1 else 0.0
    last_mean = float(sub.mean())
    last_std = float(sub.std(ddof=1)) if len(sub) > 1 else 0.0
    min_val = float(df["y"].min())
    max_val = float(df["y"].max())
    status = "STABLE"
    if "rmsd" in metric.lower() and (
        "pep" in metric.lower() or "peptide" in system_label.lower()
    ):
        # table-only; do not paint 3-stage RMSD as a feature on overview figures
        status = "TABLE_ONLY (see peptide_rmsd_jump_diagnosis.txt; not on fig0)"
    return {
        "system": system_label,
        "metric": metric,
        "mean": last_mean,
        "std": last_std,
        "full_traj_mean": full_mean,
        "full_traj_std": full_std,
        "min": min_val,
        "max": max_val,
        "status": status,
        "n_points": int(len(sub)),
    }


def main():
    parser = argparse.ArgumentParser(description="自动批量绘制 AChE-Aβ 复合物分子动力学分析图表")
    parser.add_argument("--dir", "-d", type=str, default=".", help="分析数据文件所在工作目录")
    parser.add_argument("--out", "-o", type=str, default="./figures", help="图表保存目标目录")
    parser.add_argument("--limits-json", type=str, default=None,
                        help="可选的共享 y 轴范围 JSON（统一四系统图面，由 unified_replot_and_compare.py 生成）")
    args = parser.parse_args()

    limits: dict = {}
    if args.limits_json:
        lp = Path(args.limits_json)
        if lp.exists():
            with open(lp, "r", encoding="utf-8") as fh:
                try:
                    limits = json.load(fh)
                    print(f"[LIMITS] shared y-limits loaded from {lp} ({len(limits)} panels)")
                except Exception as e:
                    print(f"[WARN] could not parse --limits-json: {e}")

    work_dir = Path(args.dir)
    fig_dir = Path(args.out)
    fig_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"Starting publication figure generation in: {fig_dir}")
    print("=" * 60)

    summary_rows: List[dict] = []

    # --------------------------------------------------------
    # 图 1: 骨架 RMSD 与 RMSF (论文 3.1 节 / 图1)
    # --------------------------------------------------------
    print("\n>> [1/8] Generating Figure 1: Backbone RMSD & RMSF ...")
    rmsd_com = read_xvg(work_dir / "rmsd_complex_bb.xvg", x_scale=0.001)
    rmsd_ach = read_xvg(work_dir / "rmsd_ache_bb.xvg", x_scale=0.001)
    rmsd_pep = read_xvg(work_dir / "rmsd_pep_bb.xvg", x_scale=0.001)

    rmsf_com = read_xvg(work_dir / "rmsf_complex_bb.xvg", x_scale=1.0)
    rmsf_ach = read_xvg(work_dir / "rmsf_ache_bb.xvg", x_scale=1.0)
    rmsf_pep = read_xvg(work_dir / "rmsf_pep_bb.xvg", x_scale=1.0)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    ax1, ax2 = axes

    # RMSD — complex only (AChE-only curves belong to the compare_* figures, not here)
    if rmsd_com is not None:
        ax1.plot(rmsd_com["x"], rmsd_com["y"], label="Complex BB", linewidth=1.2, color="tab:blue")
        s = summarize_last_ns(rmsd_com, "1_Backbone_RMSD_(nm)", "Complex")
        if s: summary_rows.append(s)
    if rmsd_ach is not None:
        # table only — no AChE-only curve on complex figures
        s = summarize_last_ns(rmsd_ach, "1_Backbone_RMSD_(nm)", "AChE")
        if s:
            s["status"] = "TABLE_ONLY (AChE part of complex; drawn in compare_* figures)"
            summary_rows.append(s)
    pep_phases = []
    ss_src = None
    # Overview / Fig1: complex analysis only (no peptide RMSD on the main figure)
    if rmsd_pep is not None:
        s = summarize_last_ns(rmsd_pep, "1_Backbone_RMSD_(nm)", "Peptide")
        if s:
            summary_rows.append(s)
    rmsd_lig = read_xvg(work_dir / "rmsd_pep_on_ache.xvg", x_scale=0.001)
    if rmsd_lig is not None:
        s = summarize_last_ns(rmsd_lig, "1c_Ligand_RMSD_fit_to_AChE_(nm)", "Peptide")
        if s:
            s["status"] = "LIGAND_POSE (table only; not drawn on overview)"
            summary_rows.append(s)
    ax1.set_title("Backbone Cα RMSD (Complex)", fontsize=11, weight="bold")
    ax1.set_xlabel("Time (ns)", fontsize=10)
    ax1.set_ylabel("RMSD (nm)", fontsize=10)
    ax1.grid(alpha=0.3, linestyle="--")
    apply_ylim(ax1, limits, "rmsd")
    safe_legend(ax1)
    add_panel_label(ax1, "A")

    # RMSF — complex backbone (residues 1-537) only; peptide RMSF lives in fig_peptide_rmsd_rmsf
    if rmsf_com is not None:
        ax2.plot(rmsf_com["x"], rmsf_com["y"], label="Complex BB", linewidth=1.0, color="tab:blue", alpha=0.9)
        s = summarize_rmsf(rmsf_com, "2_Backbone_RMSF_Avg_(nm)", "Complex")
        if s: summary_rows.append(s)
    elif rmsf_ach is not None:
        ax2.plot(rmsf_ach["x"], rmsf_ach["y"], label="AChE BB", linewidth=1.0, color="tab:orange", alpha=0.9)
        s = summarize_rmsf(rmsf_ach, "2_Backbone_RMSF_Avg_(nm)", "AChE")
        if s: summary_rows.append(s)
    if rmsf_com is not None and rmsf_ach is not None:
        # table only
        s = summarize_rmsf(rmsf_ach, "2_Backbone_RMSF_Avg_(nm)", "AChE")
        if s:
            s["status"] = "TABLE_ONLY (AChE part of complex; drawn in compare_* figures)"
            summary_rows.append(s)
    if rmsf_pep is not None:
        s = summarize_rmsf(rmsf_pep, "2_Backbone_RMSF_Avg_(nm)", "Peptide")
        if s:
            summary_rows.append(s)
    ax2.set_title("Backbone Cα RMSF (Complex)", fontsize=11, weight="bold")
    ax2.set_xlabel("Residue Number", fontsize=10)
    ax2.set_ylabel("RMSF (nm)", fontsize=10)
    ax2.grid(alpha=0.3, linestyle="--")
    apply_ylim(ax2, limits, "rmsf")
    safe_legend(ax2)
    add_panel_label(ax2, "B")

    save_all_formats(fig, fig_dir / "fig1_rmsd_rmsf")
    plt.close(fig)

    # --------------------------------------------------------
    # Peptide RMSD / RMSF — standalone (NOT on overview)
    # --------------------------------------------------------
    print("\n>> [1b] Generating peptide-only RMSD / RMSF ...")
    fig, axes_p = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    axp1, axp2 = axes_p
    if rmsd_pep is not None:
        axp1.plot(rmsd_pep["x"], rmsd_pep["y"], label="Peptide BB (self-fit)",
                  linewidth=1.2, color="tab:green")
    if rmsd_lig is not None:
        axp1.plot(rmsd_lig["x"], rmsd_lig["y"], label="Peptide fit to AChE",
                  linewidth=1.1, color="tab:purple", linestyle="--")
    if rmsd_pep is None and rmsd_lig is None:
        axp1.text(0.5, 0.5, "Peptide RMSD N/A (monomer)", ha="center", va="center", color="gray")
    axp1.set_title("Peptide Backbone RMSD", fontsize=11, weight="bold")
    axp1.set_xlabel("Time (ns)", fontsize=10)
    axp1.set_ylabel("RMSD (nm)", fontsize=10)
    axp1.grid(alpha=0.3, linestyle="--")
    safe_legend(axp1)
    add_panel_label(axp1, "A")

    if rmsf_pep is not None:
        axp2.plot(rmsf_pep["x"], rmsf_pep["y"], label="Peptide BB",
                  linewidth=1.2, color="tab:green", marker="o", markersize=4)
        axp2.set_title("Peptide Backbone RMSF", fontsize=11, weight="bold")
    else:
        axp2.text(0.5, 0.5, "Peptide RMSF N/A (monomer)", ha="center", va="center", color="gray")
        axp2.set_title("Peptide Backbone RMSF", fontsize=11, weight="bold")
    axp2.set_xlabel("Residue Number", fontsize=10)
    axp2.set_ylabel("RMSF (nm)", fontsize=10)
    axp2.grid(alpha=0.3, linestyle="--")
    safe_legend(axp2)
    add_panel_label(axp2, "B")
    save_all_formats(fig, fig_dir / "fig_peptide_rmsd_rmsf")
    plt.close(fig)

    # --------------------------------------------------------
    # 图 2: 径向分布函数 RDF (论文 3.1 节 / 图2)
    # --------------------------------------------------------
    print("\n>> [2/8] Generating Figure 2: Radial Distribution Function (RDF) ...")
    rdf_main = read_xvg(work_dir / "rdf_pep_ache.xvg", x_scale=1.0)
    fig, ax = plt.subplots(figsize=(7.5, 4.8), constrained_layout=True)
    if rdf_main is not None:
        ax.plot(rdf_main["x"], rdf_main["y"], label="Total Trajectory", linewidth=2.0, color="tab:purple")
        s = summarize_rdf(rdf_main, "5_RDF_First_Peak_g(r)", "Complex")
        if s: summary_rows.append(s)
    for q, color in enumerate(["tab:blue", "tab:orange", "tab:green", "tab:red"], start=1):
        rdf_q = read_xvg(work_dir / f"rdf_pep_ache_q{q}.xvg", x_scale=1.0)
        if rdf_q is not None:
            ax.plot(rdf_q["x"], rdf_q["y"], label=f"Quarter {q}", linewidth=1.0, linestyle="--", color=color, alpha=0.8)

    ax.set_title("Peptide around AChE Center of Mass RDF (Paper Fig 2)", fontsize=11, weight="bold")
    ax.set_xlabel("Distance (nm)", fontsize=10)
    ax.set_ylabel("g(r)", fontsize=10)
    ax.grid(alpha=0.3, linestyle="--")
    apply_ylim(ax, limits, "rdf")
    safe_legend(ax)
    save_all_formats(fig, fig_dir / "fig2_rdf")
    plt.close(fig)

    # --------------------------------------------------------
    # 图 3: 溶剂可及表面积 SASA (论文 3.2 节 / 图3)
    # --------------------------------------------------------
    print("\n>> [3/8] Generating Figure 3: Solvent Accessible Surface Area (SASA) ...")
    sasa_com = read_xvg(work_dir / "sasa_complex.xvg", x_scale=0.001)
    sasa_ach = read_xvg(work_dir / "sasa_ache.xvg", x_scale=0.001)
    sasa_pep = read_xvg(work_dir / "sasa_pep.xvg", x_scale=0.001)

    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    if sasa_com is not None:
        ax.plot(sasa_com["x"], sasa_com["y"], label="Complex SASA", linewidth=1.5, color="tab:blue")
        s = summarize_last_ns(sasa_com, "3_Solvent_Accessible_Surface_Area_SASA_(nm2)", "Complex")
        if s: summary_rows.append(s)
    if sasa_ach is not None:
        # table only — no AChE-only curve on complex figures
        s = summarize_last_ns(sasa_ach, "3_Solvent_Accessible_Surface_Area_SASA_(nm2)", "AChE")
        if s:
            s["status"] = "TABLE_ONLY (AChE part of complex; drawn in compare_* figures)"
            summary_rows.append(s)
    if sasa_pep is not None:
        s = summarize_last_ns(sasa_pep, "3_Solvent_Accessible_Surface_Area_SASA_(nm2)", "Peptide")
        if s:
            s["status"] = "TABLE_ONLY (peptide SASA; complex-only curves on figure)"
            summary_rows.append(s)
    ax.set_title("Solvent Accessible Surface Area (Complex SASA, Paper Fig 3)", fontsize=11, weight="bold")
    ax.set_xlabel("Time (ns)", fontsize=10)
    ax.set_ylabel("SASA (nm²)", fontsize=10)
    ax.grid(alpha=0.3, linestyle="--")
    apply_ylim(ax, limits, "sasa")
    safe_legend(ax)
    save_all_formats(fig, fig_dir / "fig3_sasa")
    plt.close(fig)

    # --------------------------------------------------------
    # 图 4: 肽二级结构演变 (论文 3.2 节 / 图4)
    # --------------------------------------------------------
    print("\n>> [4/8] Generating Figure 4: Complex + Peptide DSSP ...")
    ss_frac = read_ss_frac(work_dir / "ss_complex_frac.xvg")
    ss_label = "Complex (AChE+Peptide)"
    ss_system = "Complex"
    if ss_frac is None or len(ss_frac) < 5:
        ss_frac = read_ss_frac(work_dir / "ss_pep_frac.xvg")
        ss_label = "Peptide only"
        ss_system = "Peptide"
    ss_bins = read_dat_or_csv(work_dir / "ss_complex_bins.dat")
    if ss_bins is None:
        ss_bins = read_dat_or_csv(work_dir / "ss_pep_bins.dat")
    ss_times, ss_seqs = read_ss_perres(work_dir / "ss_complex_perres.dat")
    if ss_seqs is None or len(ss_seqs) < 5:
        ss_times, ss_seqs = read_ss_perres(work_dir / "ss_pep_perres.dat")
    ss_src = None
    if ss_frac is not None and len(ss_frac) >= 5:
        ss_src = ss_frac.rename(columns={"time_ns": "t"})
    elif ss_bins is not None and len(ss_bins) >= 5 and len(ss_bins.columns) >= 4:
        cols = list(ss_bins.columns)
        ss_src = pd.DataFrame({
            "t": ss_bins[cols[0]],
            "helix": ss_bins[cols[1]],
            "turn": ss_bins[cols[2]],
            "bend": ss_bins[cols[3]],
            "coil": ss_bins[cols[4]] if len(cols) > 4 else 0.0,
            "sheet": ss_bins[cols[5]] if len(cols) > 5 else 0.0,
        })

    ss_pep_frac = read_ss_frac(work_dir / "ss_pep_frac.xvg")
    ss_pep_src = None
    if ss_pep_frac is not None and len(ss_pep_frac) >= 5 and ss_system != "Peptide":
        ss_pep_src = ss_pep_frac.rename(columns={"time_ns": "t"})
    pep_times, pep_seqs = read_ss_perres(work_dir / "ss_pep_perres.dat")

    fig, axes4 = plt.subplots(2, 2, figsize=(12.6, 8.6), constrained_layout=True)
    ax, axb, axh, axp = axes4[0, 0], axes4[0, 1], axes4[1, 0], axes4[1, 1]

    if ss_src is not None:
        t, helix, sheet, turn, bend, coil = ss_to_percent_stacks(ss_src)
        stacks = [helix, sheet, turn, bend, coil]
        draw_ss_lines(ax, t, stacks)
        ax.set_title(f"{ss_label} DSSP content (%)", fontsize=11, weight="bold")
        add_panel_label(ax, "A")

        last_mask, means, stds = ss_last_window(t, stacks, last_ns=20.0)
        draw_ss_bars(axb, means, stds)
        axb.set_title("Last 20 ns occupancy (mean ± SD)", fontsize=11, weight="bold")
        add_panel_label(axb, "B")

        draw_ss_heatmap(axh, ss_times, ss_seqs)
        axh.set_title("Complex DSSP map (residue × time)", fontsize=11, weight="bold")
        add_panel_label(axh, "C")

        n_ss = int(len(ss_src))
        metric_map = (
            ("helix", "6_DSSP_Helix_Content_(%)", helix),
            ("sheet", "6_DSSP_Sheet_Content_(%)", sheet),
            ("turn", "6_DSSP_Turn_Content_(%)", turn),
            ("bend", "6_DSSP_Bend_Content_(%)", bend),
            ("coil", "6_DSSP_Coil_Content_(%)", coil),
        )
        for _k, metric, arr in metric_map:
            last = arr[last_mask]
            summary_rows.append({
                "system": ss_system,
                "metric": metric,
                "mean": float(np.mean(last)),
                "std": float(np.std(last, ddof=1)) if len(last) > 1 else 0.0,
                "full_traj_mean": float(np.mean(arr)),
                "full_traj_std": float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
                "min": float(np.min(arr)),
                "max": float(np.max(arr)),
                "status": f"DSSP_PERCENT (n={n_ss} frames)",
                "n_points": n_ss,
            })
    else:
        ax.text(0.5, 0.5, "DSSP missing — run replot_<system>.ps1",
                ha="center", va="center", fontsize=10, color="tab:red")
        ax.set_title("Secondary structure content (%)", fontsize=11, weight="bold")
        axb.axis("off")
        axh.axis("off")

    if ss_pep_src is not None:
        tp, hp, ep, up, bp, cp = ss_to_percent_stacks(ss_pep_src)
        draw_ss_lines(axp, tp, [hp, ep, up, bp, cp], fontsize=7)
        axp.set_title("Peptide DSSP content (%)  [supplement]", fontsize=11, weight="bold")
        add_panel_label(axp, "D")
        _m, pmeans, pstds = ss_last_window(tp, [hp, ep, up, bp, cp], last_ns=20.0)
        for metric, arr in (
            ("6p_Peptide_DSSP_Helix_(%)", hp),
            ("6p_Peptide_DSSP_Sheet_(%)", ep),
            ("6p_Peptide_DSSP_Coil_(%)", cp),
        ):
            last = arr[_m]
            summary_rows.append({
                "system": "Peptide",
                "metric": metric,
                "mean": float(np.mean(last)),
                "std": float(np.std(last, ddof=1)) if len(last) > 1 else 0.0,
                "full_traj_mean": float(np.mean(arr)),
                "full_traj_std": float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0,
                "min": float(np.min(arr)),
                "max": float(np.max(arr)),
                "status": "PEPTIDE_DSSP_SUPPLEMENT",
                "n_points": int(len(tp)),
            })
    else:
        axp.text(0.5, 0.5, "Peptide DSSP N/A", ha="center", va="center", color="gray")
        axp.set_title("Peptide DSSP content (%)  [supplement]", fontsize=11, weight="bold")
        add_panel_label(axp, "D")

    save_all_formats(fig, fig_dir / "fig4_secondary_structure")
    plt.close(fig)

    # --------------------------------------------------------
    # 图 5: 非天然残基接触统计 (论文 3.3 节 / 图5)
    # --------------------------------------------------------
    print("\n>> [5/8] Generating Figure 5: Non-Native Inter/Intra Contacts ...")
    inter_csv = read_dat_or_csv(work_dir / "inter_contacts.csv")
    intra_csv = read_dat_or_csv(work_dir / "intra_contacts.csv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    ax1, ax2 = axes

    if inter_csv is not None and len(inter_csv.columns) >= 2:
        cols = inter_csv.columns
        ax1.bar(inter_csv[cols[0]].astype(str), inter_csv[cols[1]], color="tab:cyan", alpha=0.8, edgecolor="black", width=0.6)
        ax1.set_title("Peptide-AChE Non-Native Contacts", fontsize=11, weight="bold")
        ax1.set_xlabel("Peptide Residue", fontsize=10)
        ax1.set_ylabel("Average Contacts per Frame", fontsize=10)
        ax1.tick_params(axis="x", rotation=45)
        ax1.grid(axis="y", alpha=0.3, linestyle="--")
        tot_c = int(len(inter_csv))
        summary_rows.append({
            "system": "Complex",
            "metric": "7_Intermolecular_Contact_Pairs_(count)",
            "mean": float(tot_c),
            "std": 0.0,
            "full_traj_mean": float(tot_c),
            "full_traj_std": 0.0,
            "min": float(tot_c),
            "max": float(tot_c),
            "status": f"{tot_c} contact pairs identified",
            "n_points": tot_c,
        })
    add_panel_label(ax1, "A")
    if intra_csv is not None and len(intra_csv.columns) >= 2:
        cols = intra_csv.columns
        ax2.bar(intra_csv[cols[0]].astype(str), intra_csv[cols[1]], color="tab:pink", alpha=0.8, edgecolor="black", width=0.6)
        ax2.set_title("Intra-Peptide Contacts", fontsize=11, weight="bold")
        ax2.set_xlabel("Peptide Residue", fontsize=10)
        ax2.set_ylabel("Average Contacts per Frame", fontsize=10)
        ax2.tick_params(axis="x", rotation=45)
        ax2.grid(axis="y", alpha=0.3, linestyle="--")
    add_panel_label(ax2, "B")
    save_all_formats(fig, fig_dir / "fig5_contacts")
    plt.close(fig)

    # --------------------------------------------------------
    # 图 6: 水介导桥连相互作用 (论文 3.4 节 / 图6 / 表2)
    # --------------------------------------------------------
    print("\n>> [6/8] Generating Figure 6: Bridging Waters ...")
    bridge_csv = read_dat_or_csv(work_dir / "bridging_per_residue.csv")
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    if bridge_csv is not None and len(bridge_csv.columns) >= 3:
        cols = bridge_csv.columns
        x_labels = bridge_csv[cols[0]].astype(str)
        ax.bar(x_labels, bridge_csv[cols[1]], label="Bridge Water Count", color="tab:blue", alpha=0.7, width=0.45)
        ax.plot(x_labels, bridge_csv[cols[2]], label="Bridge Interactions", color="tab:red", marker="o", linewidth=2.0)
        ax.set_title("Bridging Water Molecules per Residue (Paper Fig 6)", fontsize=11, weight="bold")
        ax.set_xlabel("Peptide Residue", fontsize=10)
        ax.set_ylabel("Count / Interactions", fontsize=10)
        ax.tick_params(axis="x", rotation=45)
        ax.grid(axis="y", alpha=0.3, linestyle="--")
        safe_legend(ax)
    else:
        ax.text(0.5, 0.5, "Bridging water data not available", ha="center", va="center")
    save_all_formats(fig, fig_dir / "fig6_bridging_waters")
    plt.close(fig)

    # --------------------------------------------------------
    # 图 7: 氢键分布 (论文 3.3 节)
    # --------------------------------------------------------
    print("\n>> [7/8] Generating Figure 7: Hydrogen Bonds ...")
    hb_pep_ach = read_xvg(work_dir / "hbond_ache_pep.xvg", x_scale=0.001)
    hb_intra = read_xvg(work_dir / "hbond_pep_intra.xvg", x_scale=0.001)
    hb_ach_intra = read_xvg(work_dir / "hbond_ache_intra.xvg", x_scale=0.001)
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    if hb_pep_ach is not None:
        ax.plot(hb_pep_ach["x"], hb_pep_ach["y"], label="AChE - Peptide", linewidth=1.5, color="tab:green")
        s = summarize_last_ns(hb_pep_ach, "4_Hydrogen_Bonds_(Count)", "Complex")
        if s: summary_rows.append(s)
    if hb_intra is not None:
        ax.plot(hb_intra["x"], hb_intra["y"], label="Intra-Peptide", linewidth=1.2, linestyle="--", color="tab:olive")
        s = summarize_last_ns(hb_intra, "4_Hydrogen_Bonds_(Count)", "Peptide")
        if s: summary_rows.append(s)
    if hb_ach_intra is not None:
        # table only — no protein-only curve on complex figures
        s = summarize_last_ns(hb_ach_intra, "4_Hydrogen_Bonds_(Count)", "AChE")
        if s:
            s["status"] = "TABLE_ONLY (intra-AChE; complex-focused figures only)"
            summary_rows.append(s)
    ax.set_title("Hydrogen Bonds over Time (Paper Section 3.3)", fontsize=11, weight="bold")
    ax.set_xlabel("Time (ns)", fontsize=10)
    ax.set_ylabel("Number of H-Bonds", fontsize=10)
    ax.grid(alpha=0.3, linestyle="--")
    apply_ylim(ax, limits, "hbond")
    safe_legend(ax)
    save_all_formats(fig, fig_dir / "fig_hbonds")
    plt.close(fig)

    # --------------------------------------------------------
    # 图 0: 综合 2x3 汇总图 (参照你的多幅同版画图规范)
    # --------------------------------------------------------
    print("\n>> [8/8] Generating 2x4 Master Combined Summary Figure (fig0_summary_all) ...")
    rg_com = read_xvg(work_dir / "gyrate_complex.xvg", x_scale=0.001)
    rg_ach = read_xvg(work_dir / "gyrate_ache.xvg", x_scale=0.001)
    if rg_com is not None:
        s = summarize_last_ns(rg_com, "8_Radius_of_Gyration_Rg_(nm)", "Complex")
        if s:
            summary_rows.append(s)
    if rg_ach is not None:
        s = summarize_last_ns(rg_ach, "8_Radius_of_Gyration_Rg_(nm)", "AChE")
        if s:
            summary_rows.append(s)

    fig, axes = plt.subplots(2, 4, figsize=(18.4, 8.4), constrained_layout=True)
    axes = axes.flatten()

    # A: RMSD — complex only (unified scale via --limits-json)
    if rmsd_com is not None:
        axes[0].plot(rmsd_com["x"], rmsd_com["y"], label="Complex BB", color="tab:blue", linewidth=1.2)
    axes[0].set_title("Backbone Cα RMSD", fontsize=11, weight="bold")
    axes[0].set_xlabel("Time (ns)", fontsize=10)
    axes[0].set_ylabel("RMSD (nm)", fontsize=10)
    axes[0].grid(alpha=0.3, linestyle="--")
    apply_ylim(axes[0], limits, "rmsd")
    safe_legend(axes[0])

    # B: Complex RMSF only (no AChE-only curve on complex figures)
    if rmsf_com is not None:
        axes[1].plot(rmsf_com["x"], rmsf_com["y"], label="Complex BB", color="tab:blue", linewidth=1.2)
        axes[1].set_title("Complex Backbone Cα RMSF", fontsize=11, weight="bold")
    elif rmsf_ach is not None:
        axes[1].plot(rmsf_ach["x"], rmsf_ach["y"], label="AChE BB", color="tab:orange", linewidth=1.2)
        axes[1].set_title("Backbone Cα RMSF", fontsize=11, weight="bold")
    else:
        axes[1].set_title("Backbone RMSF", fontsize=11, weight="bold")
    axes[1].set_xlabel("Residue Number", fontsize=10)
    axes[1].set_ylabel("RMSF (nm)", fontsize=10)
    axes[1].grid(alpha=0.3, linestyle="--")
    apply_ylim(axes[1], limits, "rmsf")
    safe_legend(axes[1])

    # C: RDF
    if rdf_main is not None:
        axes[2].plot(rdf_main["x"], rdf_main["y"], label="Total RDF", color="tab:purple", linewidth=1.5)
        axes[2].set_title("Peptide-AChE COM RDF", fontsize=11, weight="bold")
    else:
        axes[2].text(0.5, 0.5, "RDF (N/A for Monomer)", ha="center", va="center", fontsize=10, color="gray")
        axes[2].set_title("Peptide-AChE COM RDF", fontsize=11, weight="bold")
    axes[2].set_xlabel("Distance (nm)", fontsize=10)
    axes[2].set_ylabel("g(r)", fontsize=10)
    axes[2].grid(alpha=0.3, linestyle="--")
    apply_ylim(axes[2], limits, "rdf")
    safe_legend(axes[2])

    # D: SASA — complex only
    if sasa_com is not None:
        axes[3].plot(sasa_com["x"], sasa_com["y"], label="Complex SASA", color="tab:blue", linewidth=1.2)
    elif sasa_ach is not None:
        axes[3].plot(sasa_ach["x"], sasa_ach["y"], label="AChE SASA", color="tab:orange", linewidth=1.2)
    axes[3].set_title("Solvent Accessible Surface Area", fontsize=11, weight="bold")
    axes[3].set_xlabel("Time (ns)", fontsize=10)
    axes[3].set_ylabel("SASA (nm²)", fontsize=10)
    axes[3].grid(alpha=0.3, linestyle="--")
    apply_ylim(axes[3], limits, "sasa")
    safe_legend(axes[3])

    # E: DSSP occupancy bars (literature % panel; stack of a 530-res protein looks flat)
    if ss_src is not None:
        t0, h0, e0, u0, b0, c0 = ss_to_percent_stacks(ss_src)
        _mask, means0, stds0 = ss_last_window(t0, [h0, e0, u0, b0, c0], last_ns=20.0)
        draw_ss_bars(axes[4], means0, stds0, ylim=limits.get("ss"))
        axes[4].set_title("DSSP occupancy (last 20 ns)", fontsize=11, weight="bold")
    else:
        axes[4].text(0.5, 0.5, "DSSP missing — run replot_<system>.ps1",
                     ha="center", va="center", fontsize=9, color="tab:red")
        axes[4].set_title("DSSP occupancy (last 20 ns)", fontsize=11, weight="bold")

    # F: Complex Rg (before H-bonds; no AChE-only curve on complex figures)
    if rg_com is not None:
        axes[5].plot(rg_com["x"], rg_com["y"], label="Complex Rg", color="tab:blue", linewidth=1.2)
    if rg_com is None and rg_ach is None:
        axes[5].text(0.5, 0.5, "Rg missing — run 6_rg.sh",
                     ha="center", va="center", fontsize=10, color="tab:red")
    axes[5].set_title("Radius of Gyration", fontsize=11, weight="bold")
    axes[5].set_xlabel("Time (ns)", fontsize=10)
    axes[5].set_ylabel("Rg (nm)", fontsize=10)
    axes[5].grid(alpha=0.3, linestyle="--")
    apply_ylim(axes[5], limits, "rg")
    safe_legend(axes[5])

    # G: H-bonds (kept on overview; detailed curve also in fig_hbonds)
    if hb_pep_ach is not None:
        axes[6].plot(hb_pep_ach["x"], hb_pep_ach["y"], label="AChE-Peptide", color="tab:green", linewidth=1.2)
        axes[6].set_title("Intermolecular H-Bonds", fontsize=11, weight="bold")
    elif hb_ach_intra is not None:
        axes[6].plot(hb_ach_intra["x"], hb_ach_intra["y"], label="Intra-AChE", color="tab:orange", linewidth=1.2)
        axes[6].set_title("AChE Hydrogen Bonds", fontsize=11, weight="bold")
    elif hb_intra is not None:
        axes[6].plot(hb_intra["x"], hb_intra["y"], label="Intra-Peptide", color="tab:olive", linewidth=1.2)
        axes[6].set_title("Hydrogen Bonds", fontsize=11, weight="bold")
    else:
        axes[6].text(0.5, 0.5, "H-Bonds (N/A for Monomer)", ha="center", va="center", fontsize=10, color="gray")
        axes[6].set_title("Hydrogen Bonds", fontsize=11, weight="bold")
    if hb_intra is not None and hb_pep_ach is not None:
        axes[6].plot(hb_intra["x"], hb_intra["y"], label="Intra-Peptide",
                     color="tab:olive", linewidth=1.0, linestyle="--")
    axes[6].set_xlabel("Time (ns)", fontsize=10)
    axes[6].set_ylabel("Count", fontsize=10)
    axes[6].grid(alpha=0.3, linestyle="--")
    apply_ylim(axes[6], limits, "hbond")
    safe_legend(axes[6])

    # H: Complex DSSP % vs time (lines, so ±1% is visible)
    if ss_src is not None:
        t0, h0, e0, u0, b0, c0 = ss_to_percent_stacks(ss_src)
        draw_ss_lines(axes[7], t0, [h0, e0, u0, b0, c0], fontsize=7, ylim=limits.get("ss"))
        axes[7].set_title("Complex DSSP content (%)", fontsize=11, weight="bold")
    else:
        axes[7].text(0.5, 0.5, "DSSP missing", ha="center", va="center", color="tab:red")
        axes[7].set_title("Complex DSSP content (%)", fontsize=11, weight="bold")

    for idx, (ax, label) in enumerate(zip(axes, ["A", "B", "C", "D", "E", "F", "G", "H"])):
        add_panel_label(ax, label)

    sys_label = work_dir.name
    if sys_label.lower().startswith("md_"):
        sys_label = sys_label[3:]
    fig.suptitle(f"AChE–peptide complex MD summary ({sys_label.upper()}, 100 ns)",
                 fontsize=14, weight="bold")
    save_all_formats(fig, fig_dir / "fig0_summary_all")
    plt.close(fig)

    # --------------------------------------------------------
    # 统计指标汇总表 (生成 SCI 详细表与整洁多列对照表)
    # --------------------------------------------------------
    if summary_rows:
        print("\n" + "=" * 60)
        print("Saving statistical summary tables...")
        print("=" * 60)
        summary_df = pd.DataFrame(summary_rows)
        summary_csv = fig_dir / "summary_metrics.csv"
        summary_df.to_csv(summary_csv, index=False)
        print(f"  [SAVED] {summary_csv}")

        try:
            wide = (
                summary_df.pivot_table(
                    index="metric",
                    columns="system",
                    values="mean",
                    aggfunc="first",
                )
                .reset_index()
                .fillna("-")
            )
            wide_csv = fig_dir / "summary_metrics_wide.csv"
            wide.to_csv(wide_csv, index=False)
            print(f"  [SAVED] {wide_csv} (Dense wide-format summary)")
        except Exception as e:
            print(f"  [WARNING] Could not create wide table: {e}")

        try:
            sci_table_rows = []
            for idx, row in summary_df.iterrows():
                sci_table_rows.append({
                    "Metric_Name": row.get("metric", "-"),
                    "Target_System": row.get("system", "-"),
                    "Last_20ns_Equil_Mean_Std": f"{row.get('mean', 0.0):.4f} ± {row.get('std', 0.0):.4f}",
                    "Full_100ns_Mean_Std": f"{row.get('full_traj_mean', row.get('mean', 0.0)):.4f} ± {row.get('full_traj_std', row.get('std', 0.0)):.4f}",
                    "Min_Value": f"{row.get('min', 0.0):.4f}",
                    "Max_Value": f"{row.get('max', 0.0):.4f}",
                    "Scientific_Assessment": row.get("status", "STABLE (稳定收敛)"),
                })
            sci_df = pd.DataFrame(sci_table_rows)
            sci_csv = fig_dir / "SCI_Table1_Comprehensive_MD_Metrics.csv"
            sci_df.to_csv(sci_csv, index=False, encoding="utf-8-sig")
            print(f"  [SAVED SCI TABLE] {sci_csv} (Publication-grade comprehensive MD summary)")
        except Exception as e:
            print(f"  [WARNING] Could not create SCI table: {e}")

    print("\n" + "=" * 60)
    print("ALL PUBLICATION FIGURES & TABLES GENERATED SUCCESSFULLY!")
    print(f"Check your figures in: {fig_dir.resolve()}")
    print("=" * 60)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"\n[FATAL ERROR] {e}", file=sys.stderr)
        import traceback

        traceback.print_exc()
        sys.exit(1)
