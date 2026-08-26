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

ache 是单独 AChE 单体对照 (apo, 无肽): 它的 fig0 不含任何肽相关曲线,
AChE-肽氢键面板只画复合物系统 — ache 不参与氢键对比, 也不出现在其图例中。

用法:
    python unified_replot_and_compare.py
"""
from __future__ import annotations

import json
import re
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
from plot_all import read_ss_frac, read_xvg, ss_to_percent_stacks, xvg_has_data  # noqa: E402


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


def _last_lines(text: str, n: int = 400) -> str:
    if not text:
        return ""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    return "\n".join(lines[-n:])


def ensure_gyrate(system_dir: Path) -> bool:
    """Fill gyrate_complex.xvg for ONE system if it is missing or empty.

    v2.7.3 rewrite — the old version only checked file EXISTENCE and called
    6_rg.sh through `bash <absolute Windows path>`, which fails when the
    machine's `bash` is WSL (System32\\bash.exe) and can silently leave an
    empty/header-only xvg behind (every later run then saw the file "exists"
    and the Rg panel stayed blank). Now:

      1. a file that has no numeric data rows counts as MISSING and is removed;
      2. GROMACS is called directly (gmx.exe/gmx on PATH), groups
         Protein -> System -> Backbone -> AChE, output validated;
      3. only as a last resort falls back to 6_rg.sh via bash (path converted
         to POSIX / /mnt/<drive>/ form for WSL bash);
      4. failures print the actual gmx/bash error tail.
    """
    out = system_dir / "gyrate_complex.xvg"
    if xvg_has_data(out):
        return True
    if out.exists():
        print(f"    [Rg] {out.name} exists but has no data rows -> regenerating")
        try:
            out.unlink()
        except OSError:
            pass

    tpr = system_dir / "md.tpr"
    if not tpr.exists():
        tpr = system_dir / "md_0_1.tpr"
    traj = system_dir / "md_fit.xtc"
    if not traj.exists():
        traj = system_dir / "md.xtc"
    if not traj.exists():
        traj = system_dir / "md_0_1.xtc"
    if not tpr.exists() or not traj.exists():
        print(f"    [WARN] cannot compute Rg: missing {tpr.name} / {traj.name}")
        return False

    ndx = system_dir / "index.ndx"
    groups = ["Protein", "System", "Backbone", "AChE"] if ndx.exists() else ["Protein", "System"]

    gmx = shutil.which("gmx.exe") or shutil.which("gmx")
    if gmx:
        for group in groups:
            cmd = [gmx, "gyrate", "-s", str(tpr), "-f", str(traj),
                   "-o", str(out)]
            if ndx.exists():
                cmd += ["-n", str(ndx)]
            try:
                r = subprocess.run(
                    cmd, input=f"{group}\n", capture_output=True, text=True,
                    errors="replace", cwd=str(system_dir),
                )
            except OSError as exc:
                print(f"    [Rg] direct gmx call failed: {exc}")
                r = None
            if r is not None and xvg_has_data(out):
                print(f"    [OK] Rg filled with gmx directly (group {group})")
                return True
            tail = _last_lines((r.stderr if r is not None else "") or (r.stdout if r is not None else ""))
            print(f"    [Rg] gmx gyrate (group {group}) failed rc={r.returncode if r is not None else '?'}: {tail}")
            if out.exists():
                try:
                    out.unlink()
                except OSError:
                    pass
    else:
        print("    [Rg] gmx not found on PATH; trying 6_rg.sh via bash")

    # ---- last resort: bash 6_rg.sh (Git Bash or WSL) ----
    bash = shutil.which("bash")
    if bash:
        script = str(ANALYSIS / "6_rg.sh").replace("\\", "/")
        try:
            probe = subprocess.run([bash, "-c", "uname -r 2>/dev/null"],
                                   capture_output=True, text=True, errors="replace", timeout=30)
            uname = (probe.stdout or "").strip()
            if "microsoft" in uname.lower() or "wsl" in uname.lower():
                m = re.match(r"^([A-Za-z]):/(.*)$", script)
                if m:
                    script = f"/mnt/{m.group(1).lower()}/{m.group(2)}"
        except Exception:
            pass
        try:
            r = subprocess.run([bash, script], cwd=str(system_dir),
                               capture_output=True, text=True, errors="replace")
        except OSError as exc:
            print(f"    [Rg] bash 6_rg.sh failed to start: {exc}")
            r = None
        if r is not None and xvg_has_data(out):
            print("    [OK] Rg filled via 6_rg.sh (bash)")
            return True
        tail = _last_lines((r.stdout if r is not None else "") + "\n" + (r.stderr if r is not None else ""))
        print(f"    [WARN] 6_rg.sh failed rc={r.returncode if r is not None else '?'}: {tail}")
    else:
        print("    [WARN] no bash found on PATH -> cannot run 6_rg.sh fallback")

    print(f"    [WARN] gyrate_complex.xvg still missing for {system_dir.name}; "
          f"fig0 panel F will show 'Rg missing'")
    return False


def main():
    print("=" * 70)
    print(" UNIFIED REPLOT + COMPARE  (no MD re-run, no mdrun, no md_* deletion)")
    print("=" * 70)

    limits = compute_unified_limits()

    for s in SYSTEMS:
        d = GMX_ROOT / f"md_{s}"
        print(f"\n>> [{s}] ensuring Rg files ...")
        if not xvg_has_data(d / "gyrate_complex.xvg"):
            ok = ensure_gyrate(d)
            if not ok:
                print(f"    [WARN] {s}: gyrate_complex.xvg still missing/empty -> fig0 panel F will show a notice.")
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
        tag = "apo 对照 (无肽)" if s == "ache" else "复合物-only, 统一 y 轴"
        print(f"   gromacs_md/md_{s}/figures/fig0_summary_all.png   ({tag})")
    for _ref, c in COMPARE_PAIRS:
        print(f"   gromacs_md/compare_ache_vs_{c}/fig_compare.png + compare_summary.csv")
        print(f"      (氢键面板 F 只含 {c} 复合物的 AChE-Peptide 曲线; ache 不参与氢键对比)")
    print(f"   shared limits: {LIMITS_JSON}")
    print("=" * 70)


if __name__ == "__main__":
    main()
