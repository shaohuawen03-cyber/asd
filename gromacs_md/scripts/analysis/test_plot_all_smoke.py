#!/usr/bin/env python3
"""Reproduce the v2.6 fig4 legend crash and check the v2.6.1+ plot pack.

Also locks the v2.7.2/v2.7.3 apo-control behavior:
  - apo vs complex is decided by the pdb2gmx topology (chain B itp), which is
    what was actually simulated;
  - AChE-peptide H-bond comparison exists for complexes only: apo systems
    never get an "AChE-Peptide" curve or legend entry;
  - header-only/empty xvg files count as MISSING (xvg_has_data);
  - the Rg panel is never silently blank (notice shown when data is missing).

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
import plot_compare_systems  # noqa: E402


def _write_xvg(path: Path, xs, ys, comment="# dummy"):
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f"# {comment}\n")
        fh.write("@    title \"dummy\"\n")
        for x, y in zip(xs, ys):
            fh.write(f"{x:.4f}  {y:.6f}\n")


def _write_itp(path: Path, natoms: int):
    """Minimal pdb2gmx-style chain itp with a [ atoms ] section."""
    with path.open("w", encoding="utf-8") as fh:
        fh.write("; dummy chain itp\n")
        fh.write("[ moleculetype ]\nProtein_chain_A 3\n")
        fh.write("[ atoms ]\n; nr type resnr residue atom cgnr charge mass\n")
        for i in range(1, natoms + 1):
            fh.write(f"{i:6d} C 1 ALA CA {i} 0.0 12.01\n")


def _write_topol_files(root: Path, has_chain_b: bool, nb_atoms: int = 122):
    """pdb2gmx topology: topol.top + chain itps (chain B only for complexes)."""
    _write_itp(root / "topol_Protein_chain_A.itp", 8145)
    (root / "topol.top").write_text(
        "#include \"topol_Protein_chain_A.itp\"\n"
        + ("#include \"topol_Protein_chain_B.itp\"\n" if has_chain_b else ""),
        encoding="utf-8",
    )
    if has_chain_b:
        _write_itp(root / "topol_Protein_chain_B.itp", nb_atoms)


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
    root.mkdir(parents=True, exist_ok=True)
    t_ps = np.linspace(0.0, 100000.0, 80)  # 0-100 ns in ps
    rmsd = 0.15 + 0.02 * np.sin(np.linspace(0, 6, t_ps.size))
    rmsf_x = np.arange(1, 531)
    rmsf_y = 0.08 + 0.02 * np.sin(rmsf_x / 20.0)
    rmsf_cx = np.arange(1, 538)
    rmsf_cy = 0.08 + 0.02 * np.sin(rmsf_cx / 20.0)
    _write_xvg(root / "rmsd_complex_bb.xvg", t_ps, rmsd)
    _write_xvg(root / "rmsd_ache_bb.xvg", t_ps, rmsd * 0.95)
    _write_xvg(root / "rmsd_pep_bb.xvg", t_ps, 0.05 + 0.20 * (t_ps / t_ps.max()))
    _write_xvg(root / "rmsf_complex_bb.xvg", rmsf_cx, rmsf_cy)
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
    # complex index: contains a [ Peptide ] group (apo detection is index-driven)
    (root / "index.ndx").write_text(
        "[ AChE ]\n[ Peptide ]\n[ AChE_Backbone ]\n[ Peptide_Backbone ]\n",
        encoding="utf-8",
    )
    _write_topol_files(root, has_chain_b=True)
    return root


def _make_dummy_apo_workdir(root: Path) -> Path:
    """Standalone AChE control (ache): protein-only products, NO peptide files.

    hbond_ache_pep.xvg / rdf_pep_ache.xvg / rmsd_pep_bb.xvg etc. must be
    absent, and index.ndx must have no [ Peptide ] group.
    """
    root.mkdir(parents=True, exist_ok=True)
    t_ps = np.linspace(0.0, 100000.0, 80)  # 0-100 ns in ps
    rmsd = 0.15 + 0.02 * np.sin(np.linspace(0, 6, t_ps.size))
    rmsf_x = np.arange(1, 531)
    rmsf_y = 0.08 + 0.02 * np.sin(rmsf_x / 20.0)
    _write_xvg(root / "rmsd_complex_bb.xvg", t_ps, rmsd)   # whole protein = AChE
    _write_xvg(root / "rmsf_complex_bb.xvg", rmsf_x, rmsf_y)
    _write_xvg(root / "sasa_complex.xvg", t_ps, 214 + 3 * np.sin(np.linspace(0, 8, t_ps.size)))
    _write_xvg(root / "gyrate_complex.xvg", t_ps, 2.30 + 0.02 * np.sin(np.linspace(0, 5, t_ps.size)))
    _write_xvg(root / "hbond_ache_intra.xvg", t_ps, 380 + 6 * np.sin(np.linspace(0, 9, t_ps.size)))

    t_ns = t_ps * 0.001
    with (root / "ss_complex_frac.xvg").open("w", encoding="utf-8") as fh:
        fh.write("# time_ns helix turn bend sheet coil ppii break\n")
        for t in t_ns:
            fh.write(f"{t:.4f} 0.3428 0.1231 0.1284 0.1698 0.2091 0.0249 0.0019\n")
    with (root / "ss_complex_perres.dat").open("w", encoding="utf-8") as fh:
        fh.write("# time_ns ss_string\n")
        ss = ("H" * 18) + ("E" * 8) + ("T" * 6) + ("S" * 6) + ("~" * 12)
        for t in t_ns:
            fh.write(f"{t:.4f} {ss}\n")

    # apo index: NO [ Peptide ] group
    (root / "index.ndx").write_text(
        "[ AChE ]\n[ AChE_Backbone ]\n", encoding="utf-8"
    )
    _write_topol_files(root, has_chain_b=False)
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
        # v2.7: complex-only figures — no AChE-only overlay curves on fig0
        assert "AChE BB" not in svg0, "fig0 must not contain AChE-only curves (complex-only)"
        assert "AChE Rg" not in svg0, "fig0 must not contain AChE-only Rg curve (complex-only)"
        assert "Complex BB" in svg0
        assert "Peptide Backbone RMSD" in svg_pep
        assert "Peptide Backbone RMSF" in svg_pep
        assert "Last 20 ns occupancy" in svg4
        assert "Complex DSSP map" in svg4
        assert "Peptide DSSP content" in svg4

        sci = (out / "SCI_Table1_Comprehensive_MD_Metrics.csv").read_text(encoding="utf-8")
        assert "8_Radius_of_Gyration_Rg_(nm)" in sci
        assert "6_DSSP_Helix_Content_(%)" in sci
        assert "INDUCED_FIT_3PHASE" not in sci


def test_rmsf_profile_splits_at_chain_numbering_restart():
    """GROMACS restarts residue numbering per chain (AChE 4-542, peptide 1-7).

    Plotting that as one polyline draws a straight diagonal line from (542, y)
    back to (1, y). plot_rmsf_profile must split the profile into separate
    segments and renumber the peptide chain after the AChE chain — also when
    chain A itself contains numbering gaps (missing residues, e.g. 259-264).
    """
    # realistic complex profile: chain A 4-258, gap, 265-494, then peptide 1-7
    x = np.concatenate([
        np.arange(4, 259, dtype=float),
        np.arange(265, 495, dtype=float),
        np.arange(1, 8, dtype=float),
    ])
    y = 0.1 + 0.001 * x

    segs = plot_all.rmsf_segments(x, y)
    assert len(segs) == 3, f"gaps + chain restart must split into 3 segments, got {len(segs)}"
    assert segs[0][0][-1] == 258.0 and segs[1][0][0] == 265.0 and segs[2][0][0] == 1.0

    fig, ax = plt.subplots()
    plot_all.plot_rmsf_profile(ax, x, y, label="Complex BB", color="tab:blue")
    lines = ax.get_lines()
    assert len(lines) == 3, "three separate segments, no artificial connector lines"
    xs_pep = lines[2].get_xdata()
    # xmax after chain A (with gaps) is 494 -> peptide renumbered 495..501
    assert xs_pep[0] == 495.0 and xs_pep[-1] == 501.0, \
        f"peptide chain must continue after chain A, got {xs_pep[0]}..{xs_pep[-1]}"
    # middle gap segment keeps its real residue numbers (265-494)
    assert lines[1].get_xdata()[0] == 265.0 and lines[1].get_xdata()[-1] == 494.0
    assert len(ax.get_legend_handles_labels()[0]) == 1, "one legend entry for all segments"
    plt.close(fig)


def test_plot_all_apo_control_no_ache_peptide_hbond():
    """ache = standalone AChE control: no AChE-Peptide H-bond curve anywhere,
    labels say 'AChE BB' (not 'Complex BB'), suptitle says 'apo control',
    summary tables use system label 'AChE' (not 'Complex')."""
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        out = work / "figures"
        _make_dummy_apo_workdir(work)
        argv = ["plot_all.py", "--dir", str(work), "--out", str(out)]
        old = sys.argv
        try:
            sys.argv = argv
            plot_all.main()
        finally:
            sys.argv = old

        svg0 = (out / "fig0_summary_all.svg").read_text(encoding="utf-8", errors="ignore")
        assert "apo control" in svg0
        assert "AChE BB" in svg0
        assert "Complex BB" not in svg0
        assert "AChE SASA" in svg0
        assert "AChE Rg" in svg0
        assert "Intra-AChE" in svg0
        assert "AChE-Peptide" not in svg0, "apo control must not show AChE-Peptide H-bonds"

        svg_hb = (out / "fig_hbonds.svg").read_text(encoding="utf-8", errors="ignore")
        assert "Intra-AChE" in svg_hb
        assert "AChE - Peptide" not in svg_hb

        svg1 = (out / "fig1_rmsd_rmsf.svg").read_text(encoding="utf-8", errors="ignore")
        assert "AChE BB" in svg1
        assert "Complex BB" not in svg1

        sm = pd.read_csv(out / "summary_metrics.csv")
        assert set(sm["system"].unique()) == {"AChE"}, \
            f"apo tables must use system label 'AChE', got {set(sm['system'].unique())}"


def test_is_apo_topology_driven_and_xvg_validation():
    """is_apo must follow the pdb2gmx topology (chain B itp), and
    xvg_has_data must reject header-only/empty xvg files."""
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        complex_dir = base / "complex"
        apo_dir = base / "apo"
        _make_dummy_workdir(complex_dir)
        _make_dummy_apo_workdir(apo_dir)

        # complex: chain B itp present -> NOT apo, even if index were stale
        assert plot_all.is_apo(complex_dir) is False
        # apo: no chain B itp -> apo, even if a stale index claimed a peptide
        (apo_dir / "index.ndx").write_text(
            "[ AChE ]\n[ Peptide ]\n[ Peptide_Backbone ]\n", encoding="utf-8"
        )
        assert plot_all.is_apo(apo_dir) is True

        # header-only xvg = missing data
        empty = base / "empty.xvg"
        empty.write_text("# comment\n@ title \"x\"\n", encoding="utf-8")
        assert plot_all.xvg_has_data(empty) is False
        good = base / "good.xvg"
        _write_xvg(good, [0.0, 1.0, 2.0], [1.0, 2.0, 3.0])
        assert plot_all.xvg_has_data(good) is True


def test_plot_all_rg_panel_never_blank_when_complex_rg_missing():
    """If gyrate_complex.xvg is missing but gyrate_ache.xvg exists, fig0 panel F
    must still show the 'Rg missing' notice — never a silently blank panel."""
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        out = work / "figures"
        _make_dummy_workdir(work)   # complex dummy has both gyrate files
        (work / "gyrate_complex.xvg").unlink()
        argv = ["plot_all.py", "--dir", str(work), "--out", str(out)]
        old = sys.argv
        try:
            sys.argv = argv
            plot_all.main()
        finally:
            sys.argv = old
        svg0 = (out / "fig0_summary_all.svg").read_text(encoding="utf-8", errors="ignore")
        assert "Rg missing" in svg0, "fig0 panel F must show a notice, not be blank"
        assert "Complex Rg" not in svg0


def test_compare_apo_vs_complex_hbond_panel_is_complex_only():
    """ache (apo) vs complex: the AChE-Peptide H-bond panel F and its legend
    contain ONLY the complex curve; apo columns of the hbond CSV row are empty."""
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp)
        ref = base / "ref_apo"
        cmp_ = base / "cmp_complex"
        out = base / "compare_out"
        _make_dummy_apo_workdir(ref)
        _make_dummy_workdir(cmp_)

        argv = ["plot_compare_systems.py", "--ref", str(ref), "--cmp", str(cmp_),
                "--ref-name", "ache", "--cmp-name", "alllhrc", "--out", str(out)]
        old = sys.argv
        try:
            sys.argv = argv
            plot_compare_systems.main()
        finally:
            sys.argv = old

        svg = (out / "fig_compare.svg").read_text(encoding="utf-8", errors="ignore")
        assert "ache AChE-Peptide" not in svg, \
            "apo ref must NOT appear in the AChE-Peptide H-bond comparison"
        assert "alllhrc AChE-Peptide" in svg, "complex H-bond curve must be drawn"
        assert "ache AChE BB (apo)" in svg, "apo legend label must be 'ache AChE BB (apo)'"
        assert "ache Complex BB" not in svg
        assert "AChE apo control" in svg

        csv = pd.read_csv(out / "compare_summary.csv")
        row = csv.loc[csv["Metric"] == "AChE-Peptide_Hbonds_last20ns_(count)"].iloc[0]
        assert pd.isna(row["ache_mean"]), "apo hbond cells must be empty"
        assert pd.isna(row["Delta_alllhrc_minus_ache"]), "no hbond delta against apo"
        assert float(row["alllhrc_mean"]) > 0.0

        rdf = csv.loc[csv["Metric"] == "RDF_peak_g(r)"].iloc[0]
        assert pd.isna(rdf["ache_mean"]), "apo has no RDF"

        pep = csv.loc[csv["Metric"] == "Peptide_self-fit_RMSD_last20ns_(nm)"].iloc[0]
        assert pd.isna(pep["ache_mean"]), "apo has no peptide RMSD"


def main() -> int:
    tests = [
        test_safe_legend_accepts_fontsize,
        test_ss_percent_stacks_sum_to_100,
        test_plot_all_does_not_crash_and_writes_fig0_rg_and_dssp_percent,
        test_rmsf_profile_splits_at_chain_numbering_restart,
        test_is_apo_topology_driven_and_xvg_validation,
        test_plot_all_rg_panel_never_blank_when_complex_rg_missing,
        test_plot_all_apo_control_no_ache_peptide_hbond,
        test_compare_apo_vs_complex_hbond_panel_is_complex_only,
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
