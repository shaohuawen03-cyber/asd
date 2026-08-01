#!/usr/bin/env python3
"""
自动读取 GROMACS 分子动力学模拟与分析结果，并参照论文《乙酰胆碱酯酶-β-淀粉样肽复合物的分子动力学模拟》
批量生成出版级矢量图与位图 (SVG / PNG / PDF) 以及统计指标表。

设计说明:
  - 采用可编辑文本矢量图规范 (svg.fonttype = "none", pdf.fonttype = 42)
  - 自动转换时间单位 (GROMACS xvg 默认 ps -> 缩放 0.001 -> ns)
  - 批量生成图 1 至图 6 全部单图，以及一张 2x3 综合汇总分析图 (fig0_summary_all)
  - 导出统计汇总表格 summary_metrics.csv 及宽表 summary_metrics_wide.csv

用法:
    python3 plot_all.py [--dir 工作目录, 默认 .] [--out 图表保存目录, 默认 ./figures]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
    return {
        "system": system_label,
        "metric": metric,
        "mean": float(sub.mean()),
        "std": float(sub.std(ddof=1)) if len(sub) > 1 else 0.0,
        "n_points": int(len(sub)),
    }


def main():
    parser = argparse.ArgumentParser(description="自动批量绘制 AChE-Aβ 复合物分子动力学分析图表")
    parser.add_argument("--dir", "-d", type=str, default=".", help="分析数据文件所在工作目录")
    parser.add_argument("--out", "-o", type=str, default="./figures", help="图表保存目标目录")
    args = parser.parse_args()

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

    # RMSD
    if rmsd_com is not None:
        ax1.plot(rmsd_com["x"], rmsd_com["y"], label="Complex BB", linewidth=1.2, color="tab:blue")
        s = summarize_last_ns(rmsd_com, "rmsd_complex", "Complex")
        if s: summary_rows.append(s)
    if rmsd_ach is not None:
        ax1.plot(rmsd_ach["x"], rmsd_ach["y"], label="AChE BB", linewidth=1.2, linestyle="--", color="tab:orange")
        s = summarize_last_ns(rmsd_ach, "rmsd_ache", "AChE")
        if s: summary_rows.append(s)
    if rmsd_pep is not None:
        ax1.plot(rmsd_pep["x"], rmsd_pep["y"], label="Peptide BB", linewidth=1.5, color="tab:green")
        s = summarize_last_ns(rmsd_pep, "rmsd_pep", "Peptide")
        if s: summary_rows.append(s)
    ax1.set_title("Backbone Cα RMSD", fontsize=11, weight="bold")
    ax1.set_xlabel("Time (ns)", fontsize=10)
    ax1.set_ylabel("RMSD (nm)", fontsize=10)
    ax1.grid(alpha=0.3, linestyle="--")
    ax1.legend(frameon=False, fontsize=9)
    add_panel_label(ax1, "A")

    # RMSF
    if rmsf_pep is not None:
        ax2.plot(rmsf_pep["x"], rmsf_pep["y"], label="Peptide BB", linewidth=1.5, color="tab:green", marker="o", markersize=3)
        s = summarize_last_ns(rmsf_pep, "rmsf_pep", "Peptide")
        if s: summary_rows.append(s)
    if rmsf_ach is not None:
        ax2.plot(rmsf_ach["x"], rmsf_ach["y"], label="AChE BB", linewidth=1.0, linestyle="--", color="tab:orange", alpha=0.7)
    ax2.set_title("Backbone Cα RMSF", fontsize=11, weight="bold")
    ax2.set_xlabel("Residue Number", fontsize=10)
    ax2.set_ylabel("RMSF (nm)", fontsize=10)
    ax2.grid(alpha=0.3, linestyle="--")
    ax2.legend(frameon=False, fontsize=9)
    add_panel_label(ax2, "B")

    save_all_formats(fig, fig_dir / "fig1_rmsd_rmsf")
    plt.close(fig)

    # --------------------------------------------------------
    # 图 2: 径向分布函数 RDF (论文 3.1 节 / 图2)
    # --------------------------------------------------------
    print("\n>> [2/8] Generating Figure 2: Radial Distribution Function (RDF) ...")
    rdf_main = read_xvg(work_dir / "rdf_pep_ache.xvg", x_scale=1.0)
    fig, ax = plt.subplots(figsize=(7.5, 4.8), constrained_layout=True)
    if rdf_main is not None:
        ax.plot(rdf_main["x"], rdf_main["y"], label="Total Trajectory", linewidth=2.0, color="tab:purple")
        s = summarize_last_ns(rdf_main, "rdf_peak", "Complex")
        if s: summary_rows.append(s)
    for q, color in enumerate(["tab:blue", "tab:orange", "tab:green", "tab:red"], start=1):
        rdf_q = read_xvg(work_dir / f"rdf_pep_ache_q{q}.xvg", x_scale=1.0)
        if rdf_q is not None:
            ax.plot(rdf_q["x"], rdf_q["y"], label=f"Quarter {q}", linewidth=1.0, linestyle="--", color=color, alpha=0.8)

    ax.set_title("Peptide around AChE Center of Mass RDF (Paper Fig 2)", fontsize=11, weight="bold")
    ax.set_xlabel("Distance (nm)", fontsize=10)
    ax.set_ylabel("g(r)", fontsize=10)
    ax.grid(alpha=0.3, linestyle="--")
    ax.legend(frameon=False, fontsize=9)
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
        ax.plot(sasa_com["x"], sasa_com["y"], label="Complex Total", linewidth=1.5, color="tab:blue")
        s = summarize_last_ns(sasa_com, "sasa_complex", "Complex")
        if s: summary_rows.append(s)
    if sasa_ach is not None:
        ax.plot(sasa_ach["x"], sasa_ach["y"], label="AChE", linewidth=1.2, linestyle="--", color="tab:orange")
    if sasa_pep is not None:
        ax.plot(sasa_pep["x"], sasa_pep["y"], label="Peptide", linewidth=1.5, color="tab:green")
        s = summarize_last_ns(sasa_pep, "sasa_pep", "Peptide")
        if s: summary_rows.append(s)
    ax.set_title("Solvent Accessible Surface Area (SASA, Paper Fig 3)", fontsize=11, weight="bold")
    ax.set_xlabel("Time (ns)", fontsize=10)
    ax.set_ylabel("SASA (nm²)", fontsize=10)
    ax.grid(alpha=0.3, linestyle="--")
    ax.legend(frameon=False, fontsize=9)
    save_all_formats(fig, fig_dir / "fig3_sasa")
    plt.close(fig)

    # --------------------------------------------------------
    # 图 4: 肽二级结构演变 (论文 3.2 节 / 图4)
    # --------------------------------------------------------
    print("\n>> [4/8] Generating Figure 4: Peptide Secondary Structure Evolution ...")
    ss_bins = read_dat_or_csv(work_dir / "ss_pep_bins.dat")
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    if ss_bins is not None and len(ss_bins.columns) >= 4:
        cols = ss_bins.columns
        ax.plot(ss_bins[cols[0]], ss_bins[cols[1]], label="Helix", linewidth=1.8, marker="o", color="tab:blue")
        ax.plot(ss_bins[cols[0]], ss_bins[cols[2]], label="Turn", linewidth=1.8, marker="s", color="tab:orange")
        ax.plot(ss_bins[cols[0]], ss_bins[cols[3]], label="Bend", linewidth=1.8, marker="^", color="tab:green")
        ax.set_title("Peptide Secondary Structure Fraction over Time (Paper Fig 4)", fontsize=11, weight="bold")
        ax.set_xlabel("Time Window (ns)", fontsize=10)
        ax.set_ylabel("Fraction", fontsize=10)
        ax.grid(alpha=0.3, linestyle="--")
        ax.legend(frameon=False, fontsize=9)
    else:
        ax.text(0.5, 0.5, "Secondary structure data not available", ha="center", va="center")
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
        ax.legend(frameon=False, fontsize=9)
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
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    if hb_pep_ach is not None:
        ax.plot(hb_pep_ach["x"], hb_pep_ach["y"], label="AChE - Peptide", linewidth=1.5, color="tab:green")
        s = summarize_last_ns(hb_pep_ach, "hbond_ache_pep", "Complex")
        if s: summary_rows.append(s)
    if hb_intra is not None:
        ax.plot(hb_intra["x"], hb_intra["y"], label="Intra-Peptide", linewidth=1.2, linestyle="--", color="tab:olive")
        s = summarize_last_ns(hb_intra, "hbond_intra", "Peptide")
        if s: summary_rows.append(s)
    ax.set_title("Hydrogen Bonds over Time (Paper Section 3.3)", fontsize=11, weight="bold")
    ax.set_xlabel("Time (ns)", fontsize=10)
    ax.set_ylabel("Number of H-Bonds", fontsize=10)
    ax.grid(alpha=0.3, linestyle="--")
    ax.legend(frameon=False, fontsize=9)
    save_all_formats(fig, fig_dir / "fig_hbonds")
    plt.close(fig)

    # --------------------------------------------------------
    # 图 0: 综合 2x3 汇总图 (参照你的多幅同版画图规范)
    # --------------------------------------------------------
    print("\n>> [8/8] Generating 2x3 Master Combined Summary Figure (fig0_summary_all) ...")
    fig, axes = plt.subplots(2, 3, figsize=(16, 8.5), constrained_layout=True)
    axes = axes.flatten()

    # Panel 1: RMSD
    if rmsd_com is not None:
        axes[0].plot(rmsd_com["x"], rmsd_com["y"], label="Complex", color="tab:blue", linewidth=1.2)
    if rmsd_pep is not None:
        axes[0].plot(rmsd_pep["x"], rmsd_pep["y"], label="Peptide", color="tab:green", linewidth=1.2)
    axes[0].set_title("Backbone Cα RMSD", fontsize=11, weight="bold")
    axes[0].set_xlabel("Time (ns)", fontsize=10)
    axes[0].set_ylabel("RMSD (nm)", fontsize=10)
    axes[0].grid(alpha=0.3, linestyle="--")
    axes[0].legend(frameon=False, fontsize=9)

    # Panel 2: RMSF
    if rmsf_pep is not None:
        axes[1].plot(rmsf_pep["x"], rmsf_pep["y"], label="Peptide BB", marker="o", color="tab:green", linewidth=1.2)
    axes[1].set_title("Peptide Backbone RMSF", fontsize=11, weight="bold")
    axes[1].set_xlabel("Residue", fontsize=10)
    axes[1].set_ylabel("RMSF (nm)", fontsize=10)
    axes[1].grid(alpha=0.3, linestyle="--")
    axes[1].legend(frameon=False, fontsize=9)

    # Panel 3: RDF
    if rdf_main is not None:
        axes[2].plot(rdf_main["x"], rdf_main["y"], label="Total RDF", color="tab:purple", linewidth=1.5)
    axes[2].set_title("Peptide-AChE COM RDF", fontsize=11, weight="bold")
    axes[2].set_xlabel("Distance (nm)", fontsize=10)
    axes[2].set_ylabel("g(r)", fontsize=10)
    axes[2].grid(alpha=0.3, linestyle="--")
    axes[2].legend(frameon=False, fontsize=9)

    # Panel 4: SASA
    if sasa_com is not None:
        axes[3].plot(sasa_com["x"], sasa_com["y"], label="Complex SASA", color="tab:blue", linewidth=1.2)
    if sasa_pep is not None:
        axes[3].plot(sasa_pep["x"], sasa_pep["y"], label="Peptide SASA", color="tab:green", linewidth=1.2)
    axes[3].set_title("Solvent Accessible Surface Area", fontsize=11, weight="bold")
    axes[3].set_xlabel("Time (ns)", fontsize=10)
    axes[3].set_ylabel("SASA (nm²)", fontsize=10)
    axes[3].grid(alpha=0.3, linestyle="--")
    axes[3].legend(frameon=False, fontsize=9)

    # Panel 5: Secondary Structure
    if ss_bins is not None and len(ss_bins.columns) >= 4:
        cols = ss_bins.columns
        axes[4].plot(ss_bins[cols[0]], ss_bins[cols[1]], label="Helix", marker="o", color="tab:blue")
        axes[4].plot(ss_bins[cols[0]], ss_bins[cols[2]], label="Turn", marker="s", color="tab:orange")
        axes[4].plot(ss_bins[cols[0]], ss_bins[cols[3]], label="Bend", marker="^", color="tab:green")
    axes[4].set_title("Secondary Structure Fractions", fontsize=11, weight="bold")
    axes[4].set_xlabel("Time (ns)", fontsize=10)
    axes[4].set_ylabel("Fraction", fontsize=10)
    axes[4].grid(alpha=0.3, linestyle="--")
    axes[4].legend(frameon=False, fontsize=9)

    # Panel 6: H-bonds
    if hb_pep_ach is not None:
        axes[5].plot(hb_pep_ach["x"], hb_pep_ach["y"], label="AChE-Peptide", color="tab:green", linewidth=1.2)
    axes[5].set_title("Intermolecular H-Bonds", fontsize=11, weight="bold")
    axes[5].set_xlabel("Time (ns)", fontsize=10)
    axes[5].set_ylabel("Count", fontsize=10)
    axes[5].grid(alpha=0.3, linestyle="--")
    axes[5].legend(frameon=False, fontsize=9)

    for idx, (ax, label) in enumerate(zip(axes, ["A", "B", "C", "D", "E", "F"])):
        add_panel_label(ax, label)

    fig.suptitle("AChE-Aβ Complex Molecular Dynamics Summary", fontsize=14, weight="bold")
    save_all_formats(fig, fig_dir / "fig0_summary_all")
    plt.close(fig)

    # --------------------------------------------------------
    # 统计指标汇总表 (参照你的表格导出逻辑)
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
            wide = summary_df.pivot_table(
                index="metric",
                columns="system",
                values="mean",
                aggfunc="first",
            ).reset_index()
            wide_csv = fig_dir / "summary_metrics_wide.csv"
            wide.to_csv(wide_csv, index=False)
            print(f"  [SAVED] {wide_csv}")
        except Exception as e:
            print(f"  [WARNING] Could not create wide table: {e}")

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
