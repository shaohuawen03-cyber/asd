#!/usr/bin/env python3
"""
自动对比 AChE 单体蛋白 (Protein) 与复合物 (Complex) 的动力学指标曲线，
参考论文 3.1 和 3.2 节及用户自定义多子图作图逻辑，生成 2x3 对比汇总大图 (protein_vs_complex.{svg,png,pdf})
及统计汇总表。

用法:
    python3 plot_compare_systems.py [--protein ../md_ache] [--complex ../md_alllhrc] [--out ../compare_figures]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 矢量图规范
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["font.size"] = 10
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["figure.dpi"] = 150


def read_xvg(path: Path, x_scale: float = 1.0) -> Optional[pd.DataFrame]:
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


def save_all_formats(fig: plt.Figure, out_base: Path) -> None:
    out_base.parent.mkdir(parents=True, exist_ok=True)
    for ext, kwargs in [("png", {"dpi": 300}), ("svg", {}), ("pdf", {})]:
        f = f"{out_base}.{ext}"
        try:
            fig.savefig(f, bbox_inches="tight", **kwargs)
            print(f"  [SAVED] {f}")
        except Exception as e:
            print(f"  [WARNING] Could not save {f}: {e}")


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.1, 1.05, label, transform=ax.transAxes,
            fontsize=14, weight="bold", va="top", ha="right")


def safe_legend(ax: plt.Axes, **kwargs) -> None:
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(frameon=False, fontsize=9, **kwargs)


def summarize_last_ns(df: Optional[pd.DataFrame], metric: str, system_label: str, last_ns: float = 20.0):
    if df is None or df.empty:
        return None
    total_span = df["x"].max() - df["x"].min()
    eff_last = min(last_ns, total_span * 0.5) if total_span > 0 else 0.0
    cutoff = max(df["x"].max() - eff_last, df["x"].min())
    sub = df.loc[df["x"] >= cutoff, "y"]
    if sub.empty:
        sub = df["y"]
    return {
        "system": system_label,
        "metric": metric,
        "mean": float(sub.mean()),
        "std": float(sub.std(ddof=1)) if len(sub) > 1 else 0.0,
        "n_points": int(len(sub)),
    }


def summarize_rmsf(df: Optional[pd.DataFrame], system_label: str):
    if df is None or df.empty:
        return None
    return {
        "system": system_label,
        "metric": "rmsf",
        "mean": float(df["y"].mean()),
        "std": float(df["y"].std(ddof=1)) if len(df) > 1 else 0.0,
        "n_points": int(len(df)),
    }


def plot_compare(ax, p_df, c_df, title, xlabel, ylabel, p_color="tab:blue", c_color="tab:orange"):
    if p_df is not None and not p_df.empty:
        ax.plot(p_df["x"], p_df["y"], linewidth=1.2, linestyle="--", label="Protein", color=p_color, alpha=0.8)
    if c_df is not None and not c_df.empty:
        ax.plot(c_df["x"], c_df["y"], linewidth=1.2, label="Complex", color=c_color, alpha=0.8)
    ax.set_title(title, fontsize=11, weight="bold")
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.grid(alpha=0.3, linewidth=0.6, linestyle="--")
    safe_legend(ax, loc="best")


def main():
    parser = argparse.ArgumentParser(description="对比绘制 AChE 单体蛋白与复合物指标")
    parser.add_argument("--protein", "-p", type=str, default="../md_ache", help="单体蛋白 MD 工作目录")
    parser.add_argument("--complex", "-c", type=str, default="../md_alllhrc", help="复合物 MD 工作目录")
    parser.add_argument("--out", "-o", type=str, default="../compare_figures", help="对比图保存目录")
    args = parser.parse_args()

    prot_dir = Path(args.protein)
    comp_dir = Path(args.complex)
    fig_dir = Path(args.out)
    fig_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print(f"Starting Protein vs Complex comparison figure generation in: {fig_dir}")
    print("=" * 60)

    # 优先使用 ache 骨架与复合物骨架进行准确对比
    metrics = [
        ("rmsd_ache_bb.xvg", "rmsd_complex_bb.xvg", "Backbone RMSD", "RMSD (nm)", "Time (ns)", "rmsd", 0.001),
        ("rmsf_ache_bb.xvg", "rmsf_complex_bb.xvg", "Backbone RMSF", "RMSF (nm)", "Residue", "rmsf", 1.0),
        ("sasa_ache.xvg", "sasa_complex.xvg", "SASA", "SASA (nm²)", "Time (ns)", "sasa", 0.001),
        ("gyrate.xvg", "gyrate.xvg", "Radius of Gyration", "Rg (nm)", "Time (ns)", "gyrate", 0.001),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)
    axes = axes.flatten()
    summary_rows: List[dict] = []
    labels = ["A", "B", "C", "D", "E", "F"]

    for idx, (ax, p_file, c_file, title, ylabel, xlabel, key, xs) in enumerate(zip(axes[:4], metrics)):
        p_df = read_xvg(prot_dir / p_file, x_scale=xs)
        c_df = read_xvg(comp_dir / c_file, x_scale=xs)
        plot_compare(ax, p_df, c_df, title, xlabel, ylabel)
        add_panel_label(ax, labels[idx])

        if key == "rmsf":
            s_p, s_c = summarize_rmsf(p_df, "Protein"), summarize_rmsf(c_df, "Complex")
        else:
            s_p, s_c = summarize_last_ns(p_df, key, "Protein"), summarize_last_ns(c_df, key, "Complex")
        for s in (s_p, s_c):
            if s: summary_rows.append(s)

    # Panel 5 (E): Hydrogen Bonds (Complex)
    ax_hb = axes[4]
    hbond_df = read_xvg(comp_dir / "hbond_ache_pep.xvg", x_scale=0.001)
    if hbond_df is not None and not hbond_df.empty:
        ax_hb.plot(hbond_df["x"], hbond_df["y"], linewidth=1.2, color="tab:green", label="Complex H-Bonds", alpha=0.8)
        ax_hb.set_title("Intermolecular Hydrogen Bonds (Complex)", fontsize=11, weight="bold")
        ax_hb.set_xlabel("Time (ns)", fontsize=10)
        ax_hb.set_ylabel("Count", fontsize=10)
        ax_hb.grid(alpha=0.3, linestyle="--")
        safe_legend(ax_hb)
        add_panel_label(ax_hb, "E")
        s_hb = summarize_last_ns(hbond_df, "hbond", "Complex")
        if s_hb: summary_rows.append(s_hb)
    else:
        ax_hb.set_title("Intermolecular Hydrogen Bonds", fontsize=11, weight="bold")
        ax_hb.text(0.5, 0.5, "H-Bond Data Not Available", ha="center", va="center", color="gray")
        add_panel_label(ax_hb, "E")

    # Panel 6 (F): hide or show RDF
    axes[5].set_visible(False)

    fig.suptitle("Protein vs Complex: Comparative MD Analysis", fontsize=14, weight="bold")
    save_all_formats(fig, fig_dir / "protein_vs_complex")
    plt.close(fig)

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_csv = fig_dir / "compare_summary_metrics.csv"
        summary_df.to_csv(summary_csv, index=False)
        print(f"  [SAVED] {summary_csv}")

    print("=" * 60)
    print("ALL COMPARISON FIGURES GENERATED SUCCESSFULLY!")
    print("=" * 60)


if __name__ == "__main__":
    main()
