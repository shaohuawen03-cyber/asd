#!/usr/bin/env python3
"""
AChE 单体对照 (ache, apo) 与 AChE-肽复合物系统（例如 alllhrc）的对照图。
输出：fig_compare.{png,svg,pdf} + compare_summary.csv（两系统 mean ± SD 及差值）。

重要: ache 是单独蛋白对照 (无肽)。AChE-肽氢键只属于复合物体系，
所以氢键面板 (F) 只画有肽的系统 — ache 绝不出现氢键对比曲线，
图例同样不含 ache 的 "AChE-Peptide" 条目。

用法:
    python3 plot_compare_systems.py --ref md_ache --cmp md_alllhrc \
        --ref-name ache --cmp-name alllhrc --out ../compare_ache_vs_alllhrc \
        [--limits-json unified_limits.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_all import (  # noqa: E402
    SS_STACK_COLORS,
    SS_STACK_LABELS,
    add_panel_label,
    is_apo,
    plot_rmsf_profile,
    read_ss_frac,
    read_xvg,
    safe_legend,
    save_all_formats,
    ss_last_window,
    ss_to_percent_stacks,
)

mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["font.size"] = 10
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["figure.dpi"] = 150

COL_REF = "#c0392b"   # ache  (red)
COL_CMP = "#2166ac"   # 对照系统 (blue)


def _has_peptide(system_dir: Path) -> bool:
    """True when the system is an AChE-peptide complex (has peptide products).

    ache is the standalone AChE control (apo, no peptide): False for it, so
    ache never contributes an AChE-peptide H-bond curve or legend entry.
    """
    return not is_apo(system_dir)


def _sys_label(name: str, has_pep: bool, suffix: str) -> str:
    """Legend label: complexes '<name> Complex BB' / apo '<name> AChE BB (apo)'."""
    if has_pep:
        return f"{name} Complex {suffix}" if suffix else name
    return f"{name} AChE {suffix} (apo)" if suffix else f"{name} (apo)"


def _load(path: Path, x_scale: float):
    return read_xvg(path, x_scale=x_scale)


def _last_ns(df: pd.DataFrame, last_ns: float = 20.0):
    """last-20-ns mean ± SD (adaptive for short trajectories)."""
    if df is None or df.empty:
        return None, None
    total_span = df["x"].max() - df["x"].min()
    eff = min(last_ns, total_span * 0.5) if total_span > 0 else 0.0
    cutoff = max(df["x"].max() - eff, df["x"].min())
    sub = df.loc[df["x"] >= cutoff, "y"]
    return float(sub.mean()), float(sub.std(ddof=1)) if len(sub) > 1 else 0.0


def _apply_ylim(ax, limits, key):
    pair = limits.get(key)
    if not pair:
        return
    try:
        lo, hi = float(pair[0]), float(pair[1])
    except (TypeError, ValueError, IndexError):
        return
    if hi > lo:
        ax.set_ylim(lo, hi)


def main():
    ap = argparse.ArgumentParser(description="两个复合物系统的对照图 (ache vs X)")
    ap.add_argument("--ref", required=True, help="参考系统目录 (e.g. md_ache)")
    ap.add_argument("--cmp", required=True, help="对照系统目录 (e.g. md_alllhrc)")
    ap.add_argument("--ref-name", required=True)
    ap.add_argument("--cmp-name", required=True)
    ap.add_argument("--out", required=True, help="输出目录 (e.g. ../compare_ache_vs_alllhrc)")
    ap.add_argument("--limits-json", default=None, help="共享 y 轴范围 (与四系统 fig0 一致)")
    args = ap.parse_args()

    ref, cmp_ = Path(args.ref), Path(args.cmp)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    limits = {}
    if args.limits_json and Path(args.limits_json).exists():
        with open(args.limits_json, "r", encoding="utf-8") as fh:
            limits = json.load(fh)

    # ---------- data ----------
    rmsd_r = _load(ref / "rmsd_complex_bb.xvg", 0.001)
    rmsd_c = _load(cmp_ / "rmsd_complex_bb.xvg", 0.001)
    rmsf_r = _load(ref / "rmsf_complex_bb.xvg", 1.0)
    rmsf_c = _load(cmp_ / "rmsf_complex_bb.xvg", 1.0)
    sasa_r = _load(ref / "sasa_complex.xvg", 0.001)
    sasa_c = _load(cmp_ / "sasa_complex.xvg", 0.001)
    rg_r = _load(ref / "gyrate_complex.xvg", 0.001)
    rg_c = _load(cmp_ / "gyrate_complex.xvg", 0.001)
    hb_r = _load(ref / "hbond_ache_pep.xvg", 0.001)
    hb_c = _load(cmp_ / "hbond_ache_pep.xvg", 0.001)
    pep_r = _load(ref / "rmsd_pep_bb.xvg", 0.001)
    pep_c = _load(cmp_ / "rmsd_pep_bb.xvg", 0.001)
    rdf_r = _load(ref / "rdf_pep_ache.xvg", 1.0)
    rdf_c = _load(cmp_ / "rdf_pep_ache.xvg", 1.0)

    ss_r = read_ss_frac(ref / "ss_complex_frac.xvg")
    ss_c = read_ss_frac(cmp_ / "ss_complex_frac.xvg")

    ref_has_pep = _has_peptide(ref)
    cmp_has_pep = _has_peptide(cmp_)

    fig, axes = plt.subplots(2, 3, figsize=(15.6, 8.8), constrained_layout=True)
    axa, axb, axc, axd, axe, axf = axes.flatten()

    # A: RMSD
    if rmsd_r is not None:
        axa.plot(rmsd_r["x"], rmsd_r["y"], label=_sys_label(args.ref_name, ref_has_pep, "BB"),
                 color=COL_REF, linewidth=1.3)
    if rmsd_c is not None:
        axa.plot(rmsd_c["x"], rmsd_c["y"], label=_sys_label(args.cmp_name, cmp_has_pep, "BB"),
                 color=COL_CMP, linewidth=1.3)
    axa.set_title("Backbone Cα RMSD")
    axa.set_xlabel("Time (ns)"); axa.set_ylabel("RMSD (nm)")
    axa.grid(alpha=0.3, linestyle="--")
    _apply_ylim(axa, limits, "rmsd"); safe_legend(axa); add_panel_label(axa, "A")

    # B: RMSF
    if rmsf_r is not None:
        plot_rmsf_profile(axb, rmsf_r["x"], rmsf_r["y"],
                          label=_sys_label(args.ref_name, ref_has_pep, "BB"),
                          color=COL_REF, alpha=0.9, linewidth=1.0)
    if rmsf_c is not None:
        plot_rmsf_profile(axb, rmsf_c["x"], rmsf_c["y"],
                          label=_sys_label(args.cmp_name, cmp_has_pep, "BB"),
                          color=COL_CMP, alpha=0.9, linewidth=1.0)
    axb.set_title("Backbone Cα RMSF (per residue)")
    axb.set_xlabel("Residue Number"); axb.set_ylabel("RMSF (nm)")
    axb.grid(alpha=0.3, linestyle="--")
    _apply_ylim(axb, limits, "rmsf"); safe_legend(axb); add_panel_label(axb, "B")

    # C: SASA
    if sasa_r is not None:
        axc.plot(sasa_r["x"], sasa_r["y"], label=_sys_label(args.ref_name, ref_has_pep, "SASA"),
                 color=COL_REF, linewidth=1.3)
    if sasa_c is not None:
        axc.plot(sasa_c["x"], sasa_c["y"], label=_sys_label(args.cmp_name, cmp_has_pep, "SASA"),
                 color=COL_CMP, linewidth=1.3)
    axc.set_title("Solvent Accessible Surface Area")
    axc.set_xlabel("Time (ns)"); axc.set_ylabel("SASA (nm²)")
    axc.grid(alpha=0.3, linestyle="--")
    _apply_ylim(axc, limits, "sasa"); safe_legend(axc); add_panel_label(axc, "C")

    # D: Rg
    if rg_r is not None:
        axd.plot(rg_r["x"], rg_r["y"], label=_sys_label(args.ref_name, ref_has_pep, "Rg"),
                 color=COL_REF, linewidth=1.3)
    if rg_c is not None:
        axd.plot(rg_c["x"], rg_c["y"], label=_sys_label(args.cmp_name, cmp_has_pep, "Rg"),
                 color=COL_CMP, linewidth=1.3)
    axd.set_title("Radius of Gyration")
    axd.set_xlabel("Time (ns)"); axd.set_ylabel("Rg (nm)")
    axd.grid(alpha=0.3, linestyle="--")
    _apply_ylim(axd, limits, "rg"); safe_legend(axd); add_panel_label(axd, "D")

    # E: DSSP last-20-ns grouped bars
    labels = SS_STACK_LABELS
    xpos = np.arange(len(labels))
    means_r = means_c = stds_r = stds_c = None
    if ss_r is not None:
        t0, h0, e0, u0, b0, c0 = ss_to_percent_stacks(ss_r.rename(columns={"time_ns": "t"}))
        _m, means_r, stds_r = ss_last_window(t0, [h0, e0, u0, b0, c0], last_ns=20.0)
    if ss_c is not None:
        t0, h0, e0, u0, b0, c0 = ss_to_percent_stacks(ss_c.rename(columns={"time_ns": "t"}))
        _m, means_c, stds_c = ss_last_window(t0, [h0, e0, u0, b0, c0], last_ns=20.0)
    w = 0.36
    if means_r is not None:
        axe.bar(xpos - w / 2, means_r, yerr=stds_r, width=w, color=COL_REF, alpha=0.85,
                capsize=2.5, edgecolor="black", linewidth=0.5,
                label=_sys_label(args.ref_name, ref_has_pep, ""))
    if means_c is not None:
        axe.bar(xpos + w / 2, means_c, yerr=stds_c, width=w, color=COL_CMP, alpha=0.85,
                capsize=2.5, edgecolor="black", linewidth=0.5,
                label=_sys_label(args.cmp_name, cmp_has_pep, ""))
    axe.set_xticks(xpos)
    axe.set_xticklabels(labels, rotation=15)
    axe.set_ylabel("Content (%)")
    axe.set_title("DSSP occupancy (last 20 ns)")
    axe.grid(axis="y", alpha=0.25, linestyle="--")
    _apply_ylim(axe, limits, "ss")
    safe_legend(axe); add_panel_label(axe, "E")

    # F: AChE-peptide H-bonds — COMPLEXES ONLY.
    # ache is the standalone apo control (no peptide): it must never appear
    # in this panel (no curve, no legend entry).
    if ref_has_pep and hb_r is not None:
        axf.plot(hb_r["x"], hb_r["y"], label=f"{args.ref_name} AChE-Peptide",
                 color=COL_REF, linewidth=1.3)
    if cmp_has_pep and hb_c is not None:
        axf.plot(hb_c["x"], hb_c["y"], label=f"{args.cmp_name} AChE-Peptide",
                 color=COL_CMP, linewidth=1.3)
    if not ref_has_pep and not cmp_has_pep:
        axf.text(0.5, 0.5, "N/A — no peptide in either system",
                 ha="center", va="center", fontsize=10, color="gray")
    axf.set_title("Intermolecular H-Bonds (AChE-Peptide, complexes only)")
    axf.set_xlabel("Time (ns)"); axf.set_ylabel("Count")
    axf.grid(alpha=0.3, linestyle="--")
    _apply_ylim(axf, limits, "hbond"); safe_legend(axf); add_panel_label(axf, "F")

    if not ref_has_pep and cmp_has_pep:
        suptitle = (f"{args.ref_name} (AChE apo control) vs "
                    f"{args.cmp_name} (AChE–peptide complex) — 100 ns each")
    elif ref_has_pep and not cmp_has_pep:
        suptitle = (f"{args.ref_name} (AChE–peptide complex) vs "
                    f"{args.cmp_name} (AChE apo control) — 100 ns each")
    else:
        suptitle = f"{args.ref_name} vs {args.cmp_name} — AChE–peptide complex comparison (100 ns each)"
    fig.suptitle(suptitle, fontsize=14, weight="bold")
    save_all_formats(fig, out / "fig_compare")
    plt.close(fig)

    # ---------- CSV ----------
    rows = []

    def add(metric, ra, rb, ca, cb, fmt="{:.4f}"):
        def cell(v):
            return fmt.format(v) if v is not None else ""
        rows.append({
            "Metric": metric,
            f"{args.ref_name}_mean": cell(ra),
            f"{args.ref_name}_std": cell(rb),
            f"{args.cmp_name}_mean": cell(ca),
            f"{args.cmp_name}_std": cell(cb),
            f"Delta_{args.cmp_name}_minus_{args.ref_name}":
                fmt.format(ca - ra) if (ra is not None and ca is not None) else "",
        })

    # peptide-only rows: filled only for systems that actually have a peptide
    for metric, dr, dc, pep_only in [
        ("Backbone_RMSD_last20ns_(nm)", rmsd_r, rmsd_c, False),
        ("Peptide_self-fit_RMSD_last20ns_(nm)", pep_r, pep_c, True),
        ("SASA_last20ns_(nm2)", sasa_r, sasa_c, False),
        ("Rg_last20ns_(nm)", rg_r, rg_c, False),
        ("AChE-Peptide_Hbonds_last20ns_(count)", hb_r, hb_c, True),
    ]:
        m1 = s1 = None
        if dr is not None and (not pep_only or ref_has_pep):
            m1, s1 = _last_ns(dr)
        m2 = s2 = None
        if dc is not None and (not pep_only or cmp_has_pep):
            m2, s2 = _last_ns(dc)
        add(metric, m1, s1, m2, s2)

    for metric, dr, dc in [("Backbone_RMSF_per-residue_mean_(nm)", rmsf_r, rmsf_c)]:
        if dr is not None and dc is not None:
            add(metric, float(dr["y"].mean()), float(dr["y"].std(ddof=1)),
                float(dc["y"].mean()), float(dc["y"].std(ddof=1)))
        else:
            add(metric, None, None, None, None)

    for metric, dr, dc in [("RDF_peak_g(r)", rdf_r, rdf_c)]:
        add(metric,
            float(dr["y"].max()) if (dr is not None and ref_has_pep) else None, None,
            float(dc["y"].max()) if (dc is not None and cmp_has_pep) else None, None,
            fmt="{:.1f}")

    if means_r is not None and means_c is not None:
        for i, lab in enumerate(labels):
            add(f"DSSP_{lab}_last20ns_(%)", means_r[i], stds_r[i], means_c[i], stds_c[i], fmt="{:.2f}")

    df_out = pd.DataFrame(rows)
    csv_path = out / "compare_summary.csv"
    df_out.to_csv(csv_path, index=False, encoding="utf-8-sig")
    print(f"  [SAVED] {csv_path}")

    print(f"\n[OK] comparison figure + table written to {out.resolve()}")
    print(f"     fig_compare.png/pdf/svg  +  compare_summary.csv")


if __name__ == "__main__":
    main()
