#!/usr/bin/env python3
"""
一键：统一重画四个 100 ns 系统（alllhrc / fllhttr / ylsllqr / ache）+ 生成三组对照图。

- 不重新模拟 MD，不启动 mdrun，不删除 md_*。
- 只做"重新作图"一步 + 对照图：
    1) 若某系统缺 gyrate_*.xvg（如 alllhrc），调用 6_rg.sh 补算 Rg（需要 md_fit.xtc + gmx）。
    2) 用四个系统的真实数据计算统一 y 轴范围（写入 unified_limits.json），
       保证四张 fig0_summary_all 面板一致、刻度一致、只画复合物曲线。
    3) plot_all.py 重画每个系统全部图（fig0/fig1..fig6/hbonds/肽RMSD-RMSF/统计表）。
    4) plot_compare_systems.py 生成三组对照图：
          compare_ache_vs_alllhrc / compare_ache_vs_fllhttr / compare_ache_vs_ylsllqr
       （fig_compare.png/pdf/svg + compare_summary.csv）

用法:
    python unified_replot_and_compare.py
"""
from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent                 # gromacs_md/scripts
GMX_ROOT = HERE.parent                                 # gromacs_md
ANALYSIS = HERE / "analysis"
SYSTEMS = ["alllhrc", "fllhttr", "ylsllqr", "ache"]
COMPARE_PAIRS = [("ache", s) for s in ("alllhrc", "fllhttr", "ylsllqr")]
LIMITS_JSON = HERE / "unified_limits.json"

sys.path.insert(0, str(ANALYSIS))
from plot_all import read_ss_frac, read_xvg, ss_to_percent_stacks  # noqa: E402


def load(path: Path, x_scale: float):
    return read_xvg(path, x_scale=x_scale)


def compute_unified_limits():
    rmsd_max = rmsf_max = rdf_max = hb_max = ss_max = 0.0
    sasa_min, sasa_max = 1e9, -1e9
    rg_min, rg_max = 1e9, -1e9
    for s in SYSTEMS:
        d = GMX_ROOT / f"md_{s}"
        for fname, key in (
            ("rmsd_complex_bb.xvg", "rmsd"),
            ("rmsf_complex_bb.xvg", "rmsf"),
            ("rdf_pep_ache.xvg", "rdf"),
            ("hbond_ache_pep.xvg", "hb"),
            ("hbond_pep_intra.xvg", "hb"),
            ("sasa_complex.xvg", "sasa"),
            ("gyrate_complex.xvg", "rg"),
        ):
            x_scale = 0.001 if key not in ("rmsf", "rdf") else 1.0
            df = load(d / fname, x_scale)
            if df is None or df.empty:
                continue
            y = df["y"].to_numpy(dtype=float)
            if key == "rmsd":
                rmsd_max = max(rmsd_max, float(y.max()))
            elif key == "rmsf":
                rmsf_max = max(rmsf_max, float(y.max()))
            elif key == "rdf":
                rdf_max = max(rdf_max, float(y.max()))
            elif key == "hb":
                hb_max = max(hb_max, float(y.max()))
            elif key == "sasa":
                sasa_min = min(sasa_min, float(y.min()))
                sasa_max = max(sasa_max, float(y.max()))
            elif key == "rg":
                rg_min = min(rg_min, float(y.min()))
                rg_max = max(rg_max, float(y.max()))
        ssf = read_ss_frac(d / "ss_complex_frac.xvg")
        if ssf is not None and len(ssf) >= 5:
            _t, h, e, u, b, c = ss_to_percent_stacks(ssf.rename(columns={"time_ns": "t"}))
            for arr in (h, e, u, b, c):
                ss_max = max(ss_max, float(np.max(arr)))

    limits = {
        "rmsd": [0.0, round(max(rmsd_max * 1.1, 0.3), 3)],
        "rmsf": [0.0, round(max(rmsf_max * 1.1, 0.6), 3)],
        "rdf": [0.0, round(max(rdf_max * 1.1, 50.0), 1)],
        "sasa": [round(sasa_min - 2.0, 1), round(sasa_max + 2.0, 1)],
        "rg": [round(rg_min - 0.04, 3), round(rg_max + 0.04, 3)],
        "hbond": [0.0, float(np.ceil(max(hb_max * 1.15, 8.0)))],
        "ss": [0.0, float(np.ceil(max(ss_max * 1.2, 40.0)))],
    }
    with open(LIMITS_JSON, "w", encoding="utf-8") as fh:
        json.dump(limits, fh, indent=2)
    print("[LIMITS] unified y-limits (shared across the 4 systems):")
    for k, v in limits.items():
        print(f"    {k:6s}: {v}")
    return limits


def ensure_gyrate(system_dir: Path) -> bool:
    if (system_dir / "gyrate_complex.xvg").exists():
        return True
    bash = shutil.which("bash")
    if not bash:
        print(f"    [WARN] no bash found on PATH -> cannot auto-run 6_rg.sh for {system_dir.name}")
        return False
    r = subprocess.run([bash, str(ANALYSIS / "6_rg.sh")], cwd=system_dir,
                       capture_output=True, text=True)
    if r.returncode != 0:
        tail = (r.stdout or "")[-1200:] + (r.stderr or "")[-1200:]
        print(f"    [WARN] 6_rg.sh failed for {system_dir.name}: {tail}")
        return False
    return (system_dir / "gyrate_complex.xvg").exists()


def main():
    print("=" * 70)
    print(" UNIFIED REPLOT + COMPARE  (no MD re-run, no mdrun, no md_* deletion)")
    print("=" * 70)

    limits = compute_unified_limits()

    for s in SYSTEMS:
        d = GMX_ROOT / f"md_{s}"
        print(f"\n>> [{s}] ensuring Rg files ...")
        if not (d / "gyrate_complex.xvg").exists():
            ok = ensure_gyrate(d)
            if not ok:
                print(f"    [WARN] {s}: gyrate_complex.xvg still missing -> fig0 panel F will show a notice.")
                print(f"           Fix: run replot_{s}.ps1 once (needs md_fit.xtc + gmx).")
        else:
            print(f"    [OK] gyrate present")
        print(f">> [{s}] plot_all.py (unified limits) ...")
        r = subprocess.run(
            [sys.executable, str(ANALYSIS / "plot_all.py"),
             "--dir", str(d), "--out", str(d / "figures"),
             "--limits-json", str(LIMITS_JSON)],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            print((r.stdout or "")[-3000:])
            print((r.stderr or "")[-3000:])
            print(f"[FAIL] plot_all.py exited {r.returncode} for {s}")
        else:
            tail = (r.stdout or "").strip().splitlines()
            print("\n".join(tail[-6:]) if tail else "[OK]")

    for ref_name, cmp_name in COMPARE_PAIRS:
        out = GMX_ROOT / f"compare_ache_vs_{cmp_name}"
        print(f"\n>> [compare] {ref_name} vs {cmp_name} -> {out}")
        r = subprocess.run(
            [sys.executable, str(ANALYSIS / "plot_compare_systems.py"),
             "--ref", str(GMX_ROOT / f"md_{ref_name}"),
             "--cmp", str(GMX_ROOT / f"md_{cmp_name}"),
             "--ref-name", ref_name, "--cmp-name", cmp_name,
             "--out", str(out),
             "--limits-json", str(LIMITS_JSON)],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            print((r.stdout or "")[-2000:])
            print((r.stderr or "")[-2000:])
            print(f"[FAIL] compare {ref_name} vs {cmp_name} exited {r.returncode}")
        else:
            print((r.stdout or "").strip().splitlines()[-3:])

    print("\n" + "=" * 70)
    print(" DONE.  Check:")
    for s in SYSTEMS:
        print(f"   gromacs_md/md_{s}/figures/fig0_summary_all.png   (复合物-only, 统一 y 轴)")
    for _ref, c in COMPARE_PAIRS:
        print(f"   gromacs_md/compare_ache_vs_{c}/fig_compare.png + compare_summary.csv")
    print(f"   shared limits: {LIMITS_JSON}")
    print("=" * 70)


if __name__ == "__main__":
    main()
