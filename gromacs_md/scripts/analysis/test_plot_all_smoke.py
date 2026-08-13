#!/usr/bin/env python3
"""Reproduce the v2.6 fig4 legend crash and check the v2.6.1 plot pack.

Run:
    python3 test_plot_all_smoke.py
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import plot_all  # noqa: E402


def _write_xvg(path: Path, xs, ys, comment="# dummy"):
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f"# {comment}\n")
        fh.write("@    title \"dummy\"\n")
        for x, y in zip(xs, ys):
            fh.write(f"{x:.4f}  {y:.6f}\n")


def test_safe_legend_accepts_fontsize():
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1], label="x")
    # This is the exact call that crashed v2.6:
    # TypeError: legend() got multiple values for keyword argument 'fontsize'
    plot_all.safe_legend(ax, loc="upper right", ncol=2, fontsize=8)
    plt.close(fig)


def test_ss_percent_stacks_sum_to_100():
    n = 20
    df = pd.DataFrame({
        "t": np.linspace(0, 100, n),
        "helix": np.full(n, 0.34),
        "sheet": np.full(n, 0.17),
        "turn": np.full(n, 0.12),
        "bend": np.full(n, 0.13),
        "coil": np.full(n, 0.21),
    })
    t, h, e, u, b, c = plot_all.ss_to_percent_stacks(df)
    stacked = h + e + u + b + c
    assert t.shape == (n,)
    assert np.allclose(stacked, 100.0, atol=1e-6)
    assert 33.0 < h.mean() < 35.0


def _make_dummy_workdir(root: Path) -> Path:
    t_ps = np.linspace(0.0, 100000.0, 80)  # 0-100 ns in ps
    rmsd = 0.15 + 0.02 * np.sin(np.linspace(0, 6, t_ps.size))
    rmsf_x = np.arange(1, 531)
    rmsf_y = 0.08 + 0.02 * np.sin(rmsf_x / 20.0)
    _write_xvg(root / "rmsd_complex_bb.xvg", t_ps, rmsd)
    _write_xvg(root / "rmsd_ache_bb.xvg", t_ps, rmsd * 0.95)
    _write_xvg(root / "rmsd_pep_bb.xvg", t_ps, 0.05 + 0.20 * (t_ps / t_ps.max()))
    _write_xvg(root / "rmsf_ache_bb.xvg", rmsf_x, rmsf_y)
    _write_xvg(root / "rmsf_pep_bb.xvg", np.arange(531, 538), np.linspace(0.05, 0.20, 7))
    _write_xvg(root / "rdf_pep_ache.xvg", np.linspace(0, 3, 60),
               np.exp(-((np.linspace(0, 3, 60) - 1.2) ** 2) / 0.05) * 10)
    _write_xvg(root / "sasa_complex.xvg", t_ps, 214 + 3 * np.sin(np.linspace(0, 8, t_ps.size)))
    _write_xvg(root / "sasa_ache.xvg", t_ps, 220 + 2 * np.sin(np.linspace(0, 8, t_ps.size)))
    _write_xvg(root / "hbond_ache_pep.xvg", t_ps, 4 + np.sin(np.linspace(0, 10, t_ps.size)))
    _write_xvg(root / "hbond_pep_intra.xvg", t_ps, 1 + 0.3 * np.sin(np.linspace(0, 7, t_ps.size)))
    _write_xvg(root / "gyrate_complex.xvg", t_ps, 2.30 + 0.02 * np.sin(np.linspace(0, 5, t_ps.size)))
    _write_xvg(root / "gyrate_ache.xvg", t_ps, 2.28 + 0.015 * np.sin(np.linspace(0, 5, t_ps.size)))

    t_ns = t_ps * 0.001
    with (root / "ss_complex_frac.xvg").open("w", encoding="utf-8") as fh:
        fh.write("# time_ns helix turn bend sheet coil ppii break\n")
        for t in t_ns:
            fh.write(f"{t:.4f} 0.3428 0.1231 0.1284 0.1698 0.2091 0.0249 0.0019\n")
    with (root / "ss_pep_frac.xvg").open("w", encoding="utf-8") as fh:
        fh.write("# time_ns helix turn bend sheet coil ppii break\n")
        for t in t_ns:
            fh.write(f"{t:.4f} 0.0000 0.1881 0.0000 0.0000 0.8052 0.0067 0.0000\n")
    # short fake complex DSSP strings (helix-rich head, coil tail)
    with (root / "ss_complex_perres.dat").open("w", encoding="utf-8") as fh:
        fh.write("# time_ns ss_string\n")
        ss = ("H" * 18) + ("E" * 8) + ("T" * 6) + ("S" * 6) + ("~" * 12)
        for t in t_ns:
            fh.write(f"{t:.4f} {ss}\n")

    (root / "inter_contacts.csv").write_text(
        "residue,contacts\nA1,12.1\nL2,8.4\nL3,6.2\n", encoding="utf-8"
    )
    (root / "intra_contacts.csv").write_text(
        "residue,contacts\nA1,3.1\nL2,2.0\n", encoding="utf-8"
    )
    (root / "bridging_per_residue.csv").write_text(
        "residue,n_waters,n_bridges\nA1,4,7\nL2,2,3\n", encoding="utf-8"
    )
    return root


def test_plot_all_does_not_crash_and_writes_fig0_rg_and_dssp_percent():
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        out = work / "figures"
        _make_dummy_workdir(work)
        argv = ["plot_all.py", "--dir", str(work), "--out", str(out)]
        old = sys.argv
        try:
            sys.argv = argv
            plot_all.main()
        finally:
            sys.argv = old

        needed = [
            "fig1_rmsd_rmsf.png",
            "fig_peptide_rmsd_rmsf.png",
            "fig_peptide_rmsd_rmsf.svg",
            "fig2_rdf.png",
            "fig3_sasa.png",
            "fig4_secondary_structure.png",
            "fig5_contacts.png",
            "fig6_bridging_waters.png",
            "fig_hbonds.png",
            "fig0_summary_all.png",
            "fig0_summary_all.svg",
            "fig4_secondary_structure.svg",
            "SCI_Table1_Comprehensive_MD_Metrics.csv",
        ]
        missing = [n for n in needed if not (out / n).exists()]
        assert not missing, f"missing figures: {missing}"

        svg0 = (out / "fig0_summary_all.svg").read_text(encoding="utf-8", errors="ignore")
        svg4 = (out / "fig4_secondary_structure.svg").read_text(encoding="utf-8", errors="ignore")
        svg_pep = (out / "fig_peptide_rmsd_rmsf.svg").read_text(encoding="utf-8", errors="ignore")
        assert "Radius of Gyration" in svg0
        assert "Intermolecular H-Bonds" in svg0
        assert "DSSP occupancy" in svg0
        assert "Complex DSSP content" in svg0
        assert "Secondary Structure Fractions" not in svg0
        assert "Peptide BB (self-fit)" not in svg0
        assert "Peptide Backbone RMSD" in svg_pep
        assert "Peptide Backbone RMSF" in svg_pep
        assert "Last 20 ns occupancy" in svg4
        assert "Complex DSSP map" in svg4
        assert "Peptide DSSP content" in svg4

        sci = (out / "SCI_Table1_Comprehensive_MD_Metrics.csv").read_text(encoding="utf-8")
        assert "8_Radius_of_Gyration_Rg_(nm)" in sci
        assert "6_DSSP_Helix_Content_(%)" in sci
        assert "INDUCED_FIT_3PHASE" not in sci


def main() -> int:
    tests = [
        test_safe_legend_accepts_fontsize,
        test_ss_percent_stacks_sum_to_100,
        test_plot_all_does_not_crash_and_writes_fig0_rg_and_dssp_percent,
    ]
    failed = 0
    for fn in tests:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
        except Exception as exc:
            failed += 1
            print(f"FAIL  {fn.__name__}: {exc}")
            import traceback
            traceback.print_exc()
    if failed:
        print(f"{failed}/{len(tests)} failed")
        return 1
    print(f"ALL {len(tests)} TESTS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
