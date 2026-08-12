#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_results.py
===============

Turn the .xvg/.dat files produced by ``scripts/analyze.sh`` into the figures of
the reference paper:

    Fig 1  RMSD (complex / receptor / peptide) and RMSF
    Fig 2  RDF of the peptide around the receptor, whole run + four quarters
    Fig 3  SASA, running average and per-block means
    Fig 4  secondary-structure propensity of the peptide
    Fig 5  contacts per peptide residue (with the receptor, and internal)
    Fig 6  water bridges per peptide residue

Requires matplotlib.  If several systems are given with ``--compare`` the
per-figure curves are overlaid so the peptides can be compared directly.

Usage
-----
    python3 scripts/plot_results.py -d work/alllhrc
    python3 scripts/plot_results.py --compare work/alllhrc work/fllhttr work/ylsllqr \\
            -o figures/
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError:
    print("ERROR: matplotlib is required.\n"
          "       pip install matplotlib", file=sys.stderr)
    sys.exit(1)


# --------------------------------------------------------------------------- #
#  xvg / dat readers
# --------------------------------------------------------------------------- #

def read_xvg(path: str) -> Tuple[List[List[float]], Dict[str, str]]:
    """Return (columns, metadata).  columns[0] is x, the rest are y series."""
    rows: List[List[float]] = []
    meta: Dict[str, str] = {}
    if not os.path.isfile(path):
        return [], meta
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("@"):
                parts = line.split(None, 2)
                if len(parts) >= 3:
                    key = parts[1]
                    val = parts[2].strip().strip('"')
                    if key in ("title", "xaxis", "yaxis"):
                        meta[key] = val
                continue
            if line.startswith("#"):
                continue
            f = line.split()
            if not f:
                continue
            try:
                rows.append([float(v) for v in f])
            except ValueError:
                continue
    if not rows:
        return [], meta
    ncol = min(len(r) for r in rows)
    cols = [[r[i] for r in rows] for i in range(ncol)]
    return cols, meta


def read_table(path: str) -> List[List[str]]:
    out = []
    if not os.path.isfile(path):
        return out
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            out.append(line.split())
    return out


def tidy_time_axis(ax) -> None:
    """Keep tick labels readable however short or long the run is."""
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=6))
    for lbl in ax.get_xticklabels():
        lbl.set_rotation(30)
        lbl.set_horizontalalignment("right")


def running_mean(y: List[float], w: int) -> List[float]:
    if w <= 1 or len(y) < w:
        return y
    out, acc = [], 0.0
    for i, v in enumerate(y):
        acc += v
        if i >= w:
            acc -= y[i - w]
        out.append(acc / min(i + 1, w))
    return out


# --------------------------------------------------------------------------- #
#  Figures
# --------------------------------------------------------------------------- #

def fig_rmsd_rmsf(systems: List[Tuple[str, str]], outdir: str) -> None:
    """Fig 1: RMSD (left) and RMSF (right) for complex, receptor, peptide."""
    rows = [("Complex_backbone", "complex"),
            ("Receptor_backbone", "receptor (AChE)"),
            ("Peptide_backbone", "peptide")]

    fig, axes = plt.subplots(3, 2, figsize=(12, 11))
    any_data = False

    for r, (grp, label) in enumerate(rows):
        axl, axr = axes[r]
        for name, d in systems:
            cols, meta = read_xvg(os.path.join(d, f"rmsd_{grp}.xvg"))
            if cols:
                any_data = True
                axl.plot(cols[0], [v * 10 for v in cols[1]], lw=1.0, label=name)
        axl.set_title(f"RMSD - {label}")
        axl.set_xlabel("time (ns)")
        axl.set_ylabel("backbone RMSD (Å)")
        axl.grid(alpha=.3)
        tidy_time_axis(axl)
        if len(systems) > 1:
            axl.legend(fontsize=8)

        if grp == "Complex_backbone":
            axr.axis("off")
            axr.text(.5, .5, "RMSF is shown for the\nreceptor and the peptide",
                     ha="center", va="center", fontsize=10, color="grey")
            continue

        for name, d in systems:
            cols, meta = read_xvg(os.path.join(d, f"rmsf_{grp}.xvg"))
            if cols:
                any_data = True
                axr.plot(cols[0], [v * 10 for v in cols[1]], lw=1.0, label=name)
        axr.set_title(f"RMSF - {label}")
        axr.set_xlabel("residue")
        axr.set_ylabel("backbone RMSF (Å)")
        axr.grid(alpha=.3)
        if len(systems) > 1:
            axr.legend(fontsize=8)

    if not any_data:
        plt.close(fig)
        return
    fig.suptitle("Figure 1 - backbone RMSD and RMSF", fontsize=13)
    fig.tight_layout()
    out = os.path.join(outdir, "fig1_rmsd_rmsf.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_rdf(systems: List[Tuple[str, str]], outdir: str) -> None:
    """Fig 2: RDF over the whole run (A) and in four quarters (B)."""
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(12, 4.5))
    any_data = False

    for name, d in systems:
        cols, _ = read_xvg(os.path.join(d, "rdf_peptide_around_receptor.xvg"))
        if cols:
            any_data = True
            axa.plot([x * 10 for x in cols[0]], cols[1], lw=1.2, label=name)
    axa.axhline(1.0, color="grey", ls="--", lw=.8)
    axa.set_title("A - whole trajectory")
    axa.set_xlabel("r (Å)")
    axa.set_ylabel("g(r)  peptide around receptor")
    axa.grid(alpha=.3)
    if len(systems) > 1:
        axa.legend(fontsize=8)

    # quarters: for one system show its four quarters, for many show quarter 1..4
    if len(systems) == 1:
        name, d = systems[0]
        for q in (1, 2, 3, 4):
            cols, _ = read_xvg(os.path.join(d, f"rdf_quarter{q}.xvg"))
            if cols:
                any_data = True
                axb.plot([x * 10 for x in cols[0]], cols[1], lw=1.0,
                         label=f"quarter {q}")
    else:
        for name, d in systems:
            for q in (1, 2, 3, 4):
                cols, _ = read_xvg(os.path.join(d, f"rdf_quarter{q}.xvg"))
                if cols:
                    any_data = True
                    axb.plot([x * 10 for x in cols[0]], cols[1], lw=.8,
                             alpha=.7, label=f"{name} q{q}")
    axb.axhline(1.0, color="grey", ls="--", lw=.8)
    axb.set_title("B - four equal time intervals")
    axb.set_xlabel("r (Å)")
    axb.set_ylabel("g(r)")
    axb.grid(alpha=.3)
    axb.legend(fontsize=7)

    if not any_data:
        plt.close(fig)
        return
    fig.suptitle("Figure 2 - RDF of the peptide around the receptor", fontsize=13)
    fig.tight_layout()
    out = os.path.join(outdir, "fig2_rdf.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_sasa(systems: List[Tuple[str, str]], outdir: str) -> None:
    """Fig 3: SASA time series + convergence (running mean)."""
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(12, 4.5))
    any_data = False

    for name, d in systems:
        cols, _ = read_xvg(os.path.join(d, "sasa_complex.xvg"))
        if not cols:
            continue
        any_data = True
        t, y = cols[0], cols[1]
        axa.plot(t, y, lw=.7, alpha=.55)
        w = max(1, len(y) // 20)
        axa.plot(t, running_mean(y, w), lw=1.6, label=name)
        # convergence: cumulative mean
        cum, acc = [], 0.0
        for i, v in enumerate(y, 1):
            acc += v
            cum.append(acc / i)
        axb.plot(t, cum, lw=1.4, label=name)

    axa.set_title("A - SASA of the complex")
    axa.set_xlabel("time (ns)")
    axa.set_ylabel("SASA (nm²)")
    axa.grid(alpha=.3)
    tidy_time_axis(axa)
    axa.legend(fontsize=8)
    axb.set_title("B - convergence (cumulative mean)")
    axb.set_xlabel("time (ns)")
    axb.set_ylabel("mean SASA (nm²)")
    axb.grid(alpha=.3)
    tidy_time_axis(axb)
    axb.legend(fontsize=8)

    if not any_data:
        plt.close(fig)
        return
    fig.suptitle("Figure 3 - solvent-accessible surface area", fontsize=13)
    fig.tight_layout()
    out = os.path.join(outdir, "fig3_sasa.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_dssp(systems: List[Tuple[str, str]], outdir: str) -> None:
    """Fig 4: secondary-structure propensity of the peptide."""
    fig, ax = plt.subplots(figsize=(9, 4.5))
    any_data = False

    for name, d in systems:
        cols, meta = read_xvg(os.path.join(d, "dssp_count.xvg"))
        if not cols or len(cols) < 2:
            continue
        any_data = True
        t = cols[0]
        for ci in range(1, len(cols)):
            ax.plot(t, cols[ci], lw=1.0,
                    label=f"{name} s{ci}" if len(systems) > 1 else f"structure {ci}")

    ax.set_title("Figure 4 - secondary-structure content of the peptide")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("number of residues")
    ax.grid(alpha=.3)
    tidy_time_axis(ax)
    ax.legend(fontsize=7, ncol=2)

    if not any_data:
        plt.close(fig)
        return
    fig.tight_layout()
    out = os.path.join(outdir, "fig4_dssp.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_contacts(systems: List[Tuple[str, str]], outdir: str) -> None:
    """Fig 5: contacts per peptide residue, with the receptor and internal."""
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(13, 4.5))
    any_data = False
    nsys = max(1, len(systems))
    width = 0.8 / nsys

    for si, (name, d) in enumerate(systems):
        inter = read_table(os.path.join(d, "contacts_peptide_receptor.dat"))
        if inter:
            any_data = True
            labels = [r[0] for r in inter]
            vals = [float(r[1]) for r in inter]
            xs = [i + si * width for i in range(len(vals))]
            axa.bar(xs, vals, width=width, label=name)
            if si == nsys - 1:
                axa.set_xticks([i + width * (nsys - 1) / 2
                                for i in range(len(labels))])
                axa.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)

        intra = read_table(os.path.join(d, "contacts_peptide_intra.dat"))
        if intra:
            any_data = True
            labels = [r[0] for r in intra]
            vals = [float(r[1]) for r in intra]
            xs = [i + si * width for i in range(len(vals))]
            axb.bar(xs, vals, width=width, label=name)
            if si == nsys - 1:
                axb.set_xticks([i + width * (nsys - 1) / 2
                                for i in range(len(labels))])
                axb.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)

    axa.set_title("A - contacts with the receptor")
    axa.set_ylabel("mean number of contacts (7 Å)")
    axa.grid(alpha=.3, axis="y")
    axa.legend(fontsize=8)
    axb.set_title("B - intramolecular contacts inside the peptide")
    axb.set_ylabel("mean number of contacts (7 Å)")
    axb.grid(alpha=.3, axis="y")
    axb.legend(fontsize=8)

    if not any_data:
        plt.close(fig)
        return
    fig.suptitle("Figure 5 - non-native contacts", fontsize=13)
    fig.tight_layout()
    out = os.path.join(outdir, "fig5_contacts.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


def fig_bridges(systems: List[Tuple[str, str]], outdir: str) -> None:
    """Fig 6: bridging waters and bridging interactions per peptide residue."""
    fig, ax = plt.subplots(figsize=(11, 4.5))
    any_data = False
    nsys = max(1, len(systems))
    width = 0.8 / (2 * nsys)

    for si, (name, d) in enumerate(systems):
        rows = read_table(os.path.join(d, "water_bridges.dat"))
        if not rows:
            continue
        any_data = True
        labels = [r[0] for r in rows]
        n_w = [float(r[1]) for r in rows]
        n_b = [float(r[2]) for r in rows]
        base = [i + si * 2 * width for i in range(len(rows))]
        ax.bar(base, n_b, width=width,
               label=f"{name} bridges" if nsys > 1 else "bridges")
        ax.bar([x + width for x in base], n_w, width=width,
               label=f"{name} waters" if nsys > 1 else "bridging waters")
        if si == nsys - 1:
            ax.set_xticks([i + width * (2 * nsys - 1) / 2
                           for i in range(len(labels))])
            ax.set_xticklabels(labels, rotation=60, ha="right", fontsize=7)

    ax.set_title("Figure 6 - water-mediated bridges between peptide and receptor")
    ax.set_ylabel("count")
    ax.grid(alpha=.3, axis="y")
    ax.legend(fontsize=8)

    if not any_data:
        plt.close(fig)
        return
    fig.tight_layout()
    out = os.path.join(outdir, "fig6_water_bridges.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"  wrote {out}")


# --------------------------------------------------------------------------- #

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Plot the analysis output")
    ap.add_argument("-d", "--dir", default=None,
                    help="a single system directory (its analysis/ is used)")
    ap.add_argument("--compare", nargs="+", default=None,
                    help="several system directories to overlay")
    ap.add_argument("-o", "--outdir", default=None,
                    help="where to write the PNGs "
                         "(default: the analysis dir, or ./figures for --compare)")
    args = ap.parse_args(argv)

    if not args.dir and not args.compare:
        ap.error("give -d/--dir or --compare")

    def resolve(d: str) -> Tuple[str, str]:
        """(label, analysis_dir)"""
        d = os.path.abspath(d)
        label = os.path.basename(d.rstrip("/"))
        cand = os.path.join(d, "analysis")
        return label, (cand if os.path.isdir(cand) else d)

    if args.compare:
        systems = [resolve(d) for d in args.compare]
        outdir = args.outdir or "figures"
    else:
        systems = [resolve(args.dir)]
        outdir = args.outdir or systems[0][1]

    os.makedirs(outdir, exist_ok=True)

    print("=" * 60)
    print("plot_results.py")
    for label, d in systems:
        print(f"  {label:<20s} <- {d}")
    print(f"  output -> {outdir}")
    print("=" * 60)

    fig_rmsd_rmsf(systems, outdir)
    fig_rdf(systems, outdir)
    fig_sasa(systems, outdir)
    fig_dssp(systems, outdir)
    fig_contacts(systems, outdir)
    fig_bridges(systems, outdir)

    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
