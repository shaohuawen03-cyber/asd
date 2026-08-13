#!/usr/bin/env python3
"""
Robust peptide secondary-structure analysis for short AChE-bound peptides.

Fixes the v2.3 DSSP failure mode:
  * GROMACS 2025 ``gmx dssp -o`` writes one SS string per frame with NO time
    column and uses ``~`` (loop) / ``=`` (break) instead of ``C``.
  * ``gmx dssp -num`` column order is H B E G I P S T = ~  (not classic do_dssp).
  * A 50 ns bin on a 100 ns trajectory collapses the whole DSSP record to
    2 points, which looks like a failed calculation (helix 0 ± 0, n=2).
  * ``-sel Peptide`` may silently analyse the 3-segment AChE receptor
    (4ey6 chain breaks) instead of the 7-mer. Strings longer than the
    peptide are automatically trimmed to the last N peptide residues.

Outputs (written next to the inputs):
  ss_pep_perres.dat   time_ns + DSSP string (one row per frame)
  ss_pep_frac.xvg     per-frame helix/turn/bend/sheet/coil/ppii/break fractions
  ss_pep_bins.dat     1 ns (adaptive) window averages for Paper Fig 4
  ss_pep_rama.dat     Ramachandran-region fractions (more sensitive for 7-mers)
  ss_pep_summary.txt  human-readable report
"""
from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

HELIX_CHARS = set("HGI")
SHEET_CHARS = set("EB")
TURN_CHARS = set("T")
BEND_CHARS = set("S")
PPII_CHARS = set("P")
BREAK_CHARS = set("=")
COIL_CHARS = set("~C -")

CATEGORIES = ("helix", "turn", "bend", "sheet", "coil", "ppii", "break")


def classify_char(c: str) -> str:
    if c in HELIX_CHARS:
        return "helix"
    if c in TURN_CHARS:
        return "turn"
    if c in BEND_CHARS:
        return "bend"
    if c in SHEET_CHARS:
        return "sheet"
    if c in PPII_CHARS:
        return "ppii"
    if c in BREAK_CHARS:
        return "break"
    return "coil"


def auto_window_ns(t0: float, t1: float, testing: bool = False) -> float:
    span = max(float(t1) - float(t0), 1e-9)
    if testing or span < 0.5:
        return max(span / 20.0, 1e-4)
    if span <= 150.0:
        return 1.0
    return 50.0


def _is_comment(line: str) -> bool:
    s = line.lstrip()
    return (not s) or s.startswith(("#", "@", ";"))


def infer_times(n: int, dt_ns: float = 0.0) -> np.ndarray:
    """Map frame index to ns. 5000-frame production is 100 ns (0.02 ns/frame), not 0.1."""
    n = int(n)
    if dt_ns and dt_ns > 0:
        return np.arange(n, dtype=float) * float(dt_ns)
    if n >= 200:
        return np.linspace(0.0, 100.0, n)
    return np.arange(n, dtype=float) * 0.1


def correct_time_axis(t) -> np.ndarray:
    """5001 frames labeled 0-500 ns means dt was wrongly 0.1 ns; production is 100 ns."""
    t = np.asarray(t, dtype=float)
    if t.size == 0:
        return t
    if t.size >= 2000 and float(t.max()) > 150.0:
        print(f">> [time-axis fix] {t.max():.1f} ns / {t.size} frames -> rescale to 0-100 ns")
        t = t * (100.0 / float(t.max()))
    return t


def trim_ss_to_peptide(ss: str, nres: int) -> str:
    """Only trim when we explicitly want the peptide tail of a long protein string."""
    ss = "".join(ch for ch in ss if not ch.isspace())
    if nres > 0 and len(ss) > nres + 2:
        return ss[-nres:]
    return ss


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_ss_dat(path: Path, nres: int = 7, dt_ns: float = 0.1):
    """Parse gmx dssp -o / do_dssp -ssdump / free-form SS strings."""
    if not path.exists():
        return None, []
    times: List[float] = []
    seqs: List[str] = []
    has_time = False
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if _is_comment(line):
                continue
            parts = line.strip().split()
            if not parts:
                continue
            t_val = None
            raw = ""
            try:
                t_val = float(parts[0])
                raw = "".join(parts[1:]) if len(parts) > 1 else ""
                if raw:
                    has_time = True
                else:
                    # A pure numeric token that is actually not an SS string
                    continue
            except ValueError:
                raw = "".join(parts)
            if not raw:
                continue
            seqs.append(trim_ss_to_peptide(raw, nres))
            times.append(t_val if (has_time and t_val is not None) else float(len(seqs) - 1))
    if not seqs:
        return None, []
    if has_time:
        t = np.asarray(times, dtype=float)
        # classic do_dssp dumps ns or ps; if max > 500 treat as ps
        if t.max() > 500:
            t = t * 0.001
        return t, seqs
    return infer_times(len(seqs), dt_ns), seqs


def _legend_key(label: str) -> str:
    s = label.strip().strip('"').strip("'").lower()
    s = s.replace("_", "-").replace(" ", "-")
    aliases = {
        "a-helix": "helix",
        "alpha-helix": "helix",
        "alpha": "helix",
        "5-helix": "helix5",
        "pi-helix": "helix5",
        "3-helix": "helix3",
        "3-10-helix": "helix3",
        "3_10-helix": "helix3",
        "b-sheet": "sheet",
        "beta-sheet": "sheet",
        "extended": "sheet",
        "b-bridge": "bridge",
        "beta-bridge": "bridge",
        "bridge": "bridge",
        "coil": "coil",
        "loop": "coil",
        "structure": "structure",
        "turn": "turn",
        "bend": "bend",
        "break": "break",
        "kappa-helix": "ppii",
        "ppii": "ppii",
        "pp-helix": "ppii",
        "polyproline": "ppii",
        "h": "helix",
        "g": "helix3",
        "i": "helix5",
        "e": "sheet",
        "b": "bridge",
        "t": "turn",
        "s": "bend",
        "p": "ppii",
        "~": "coil",
        "=": "break",
        "c": "coil",
    }
    return aliases.get(s, s)


def parse_ss_num_xvg(path: Path, nres: int = 7):
    """Parse gmx dssp -num / classic do_dssp -sc, using @ sN legend when present."""
    if not path.exists():
        return None, None
    legends: Dict[int, str] = {}
    rows: List[List[float]] = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            raw = line.rstrip("\n")
            if raw.lstrip().startswith("@"):
                # @ s0 legend "H"
                low = raw.lower()
                if "legend" in low:
                    try:
                        idx_token = raw.split()[1]  # s0
                        si = int(idx_token[1:]) if idx_token.startswith(("s", "S")) else None
                    except (IndexError, ValueError):
                        si = None
                    q = raw.split("legend", 1)[-1].strip().strip('"').strip("'")
                    if si is not None:
                        legends[si] = _legend_key(q)
                continue
            if _is_comment(raw):
                continue
            parts = raw.split()
            if len(parts) < 2:
                continue
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    if not rows:
        return None, None
    arr = np.asarray(rows, dtype=float)
    t = arr[:, 0]
    if t.max() > 500:
        t = t * 0.001
    t = correct_time_axis(t)

    ncols = arr.shape[1] - 1
    # Build per-category counts
    counts = {k: np.zeros(len(arr)) for k in CATEGORIES}

    def add(name: str, col: int) -> None:
        if 1 <= col <= ncols:
            if name == "helix" or name == "helix3" or name == "helix5":
                counts["helix"] += arr[:, col]
            elif name == "bridge":
                counts["sheet"] += arr[:, col]
            elif name in counts:
                counts[name] += arr[:, col]

    if legends:
        for si, name in legends.items():
            add(name, si + 1)
    else:
        # No legend: distinguish GROMACS 2025 (time + 10 SS types) vs classic do_dssp
        # 2025: H B E G I P S T = ~
        # classic: Structure Coil B-Sheet B-Bridge Bend Turn A-Helix 5-Helix 3-Helix
        if ncols >= 10:
            # Prefer 2025 layout when the last two columns look like break/loop
            # (loop count often equals nres for a short peptide).
            last = arr[:, -1]
            if np.nanmedian(last) >= 0.4 * nres or ncols >= 11:
                mapping = {
                    1: "helix", 2: "bridge", 3: "sheet", 4: "helix3",
                    5: "helix5", 6: "ppii", 7: "bend", 8: "turn",
                    9: "break", 10: "coil",
                }
            else:
                mapping = {
                    2: "coil", 3: "sheet", 4: "bridge", 5: "bend",
                    6: "turn", 7: "helix", 8: "helix5", 9: "helix3",
                }
            for col, name in mapping.items():
                add(name, col)
        elif ncols >= 7:
            mapping = {
                2: "coil", 3: "sheet", 4: "bridge", 5: "bend",
                6: "turn", 7: "helix",
            }
            for col, name in mapping.items():
                add(name, col)
        else:
            counts["coil"] = arr[:, 1]

    totals = np.zeros(len(arr))
    for k in CATEGORIES:
        totals += counts[k]
    totals[totals <= 0] = float(nres)
    fracs = {k: counts[k] / totals for k in CATEGORIES}
    return t, fracs


def seqs_to_fracs(seqs: Sequence[str]) -> Dict[str, np.ndarray]:
    n = len(seqs)
    out = {k: np.zeros(n) for k in CATEGORIES}
    for i, s in enumerate(seqs):
        if not s:
            out["coil"][i] = 1.0
            continue
        c = Counter(classify_char(ch) for ch in s)
        tot = float(len(s))
        for k in CATEGORIES:
            out[k][i] = c.get(k, 0) / tot
    return out


def bin_fracs(times: np.ndarray, fracs: Dict[str, np.ndarray], window_ns: float):
    t0, t1 = float(times[0]), float(times[-1])
    nwin = max(1, int(math.ceil((t1 - t0) / window_ns)))
    rows = []
    for i in range(nwin):
        lo = t0 + i * window_ns
        hi = lo + window_ns
        if i == nwin - 1:
            mask = (times >= lo) & (times <= t1 + 1e-12)
        else:
            mask = (times >= lo) & (times < hi)
        if not np.any(mask):
            continue
        rec = {"time_ns": lo + 0.5 * window_ns}
        for k in CATEGORIES:
            rec[k] = float(np.mean(fracs[k][mask]))
        rec["n_frames"] = int(np.sum(mask))
        rows.append(rec)
    return rows


# ---------------------------------------------------------------------------
# Trajectory fallback (Kabsch–Sander + Ramachandran)
# ---------------------------------------------------------------------------

def _virtual_h(c_prev: np.ndarray, n_xyz: np.ndarray, ca: np.ndarray) -> np.ndarray:
    v1 = n_xyz - c_prev
    v2 = n_xyz - ca
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < 1e-8 or n2 < 1e-8:
        return n_xyz.copy()
    hdir = v1 / n1 + v2 / n2
    hn = np.linalg.norm(hdir)
    if hn < 1e-8:
        return n_xyz.copy()
    return n_xyz + 1.01 * (hdir / hn)


def _hb_energy(o, c, n, h) -> float:
    def inv(a, b):
        d = float(np.linalg.norm(a - b))
        return 0.0 if d < 0.5 else 1.0 / d
    # Kabsch & Sander 1983, distances in Angstrom, E in kcal/mol
    return 0.42 * 0.20 * 332.0 * (inv(o, n) + inv(c, h) - inv(o, h) - inv(c, n))


def _dihedral(a, b, c, d) -> float:
    b1 = b - a
    b2 = c - b
    b3 = d - c
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1n = np.linalg.norm(n1)
    n2n = np.linalg.norm(n2)
    b2n = np.linalg.norm(b2)
    if n1n < 1e-8 or n2n < 1e-8 or b2n < 1e-8:
        return 0.0
    x = float(np.dot(n1, n2) / (n1n * n2n))
    x = max(-1.0, min(1.0, x))
    y = float(np.dot(np.cross(n1, n2), b2) / (n1n * n2n * b2n))
    return math.degrees(math.atan2(y, x))


def kabsch_sander_assign(bb: Dict[str, np.ndarray]) -> str:
    """Assign DSSP-like codes for one frame. bb[name] has shape (nres, 3) in Å."""
    nres = bb["N"].shape[0]
    if nres < 3:
        return "C" * nres
    hxyz = np.zeros((nres, 3))
    for i in range(nres):
        if i == 0:
            hxyz[i] = bb["N"][i]
        else:
            hxyz[i] = _virtual_h(bb["C"][i - 1], bb["N"][i], bb["CA"][i])
    # H-bond matrix: donor j (NH) bonded to acceptor i (CO)
    hb = np.zeros((nres, nres), dtype=bool)
    for i in range(nres):
        for j in range(nres):
            if abs(i - j) < 2:
                continue
            e = _hb_energy(bb["O"][i], bb["C"][i], bb["N"][j], hxyz[j])
            if e < -0.5:
                hb[i, j] = True
    is_turn = {3: np.zeros(nres, dtype=bool),
               4: np.zeros(nres, dtype=bool),
               5: np.zeros(nres, dtype=bool)}
    for nturn, arr in is_turn.items():
        for i in range(nres - nturn):
            if hb[i, i + nturn]:
                arr[i] = True
    ss = ["C"] * nres
    # helices: two consecutive n-turns → residues i+1 .. i+n
    for nturn, code in ((5, "I"), (4, "H"), (3, "G")):
        turns = is_turn[nturn]
        for i in range(nres - nturn - 1):
            if turns[i] and turns[i + 1]:
                for k in range(i + 1, i + nturn + 1):
                    if 0 <= k < nres and ss[k] == "C":
                        ss[k] = code
    # remaining isolated turns
    for nturn in (3, 4, 5):
        for i in range(nres - nturn):
            if is_turn[nturn][i]:
                for k in (i + 1, i + 2):
                    if 0 <= k < nres and ss[k] == "C":
                        ss[k] = "T"
    # bend: CA curvature > 70 deg
    for i in range(2, nres - 2):
        if ss[i] != "C":
            continue
        v1 = bb["CA"][i] - bb["CA"][i - 2]
        v2 = bb["CA"][i + 2] - bb["CA"][i]
        n1 = np.linalg.norm(v1)
        n2 = np.linalg.norm(v2)
        if n1 < 1e-6 or n2 < 1e-6:
            continue
        ang = math.degrees(math.acos(max(-1.0, min(1.0, float(np.dot(v1, v2) / (n1 * n2))))))
        if ang > 70.0:
            ss[i] = "S"
    return "".join(ss)


def rama_assign(bb: Dict[str, np.ndarray]) -> str:
    nres = bb["N"].shape[0]
    ss = ["C"] * nres
    for i in range(1, nres - 1):
        phi = _dihedral(bb["C"][i - 1], bb["N"][i], bb["CA"][i], bb["C"][i])
        psi = _dihedral(bb["N"][i], bb["CA"][i], bb["C"][i], bb["N"][i + 1])
        if -90.0 <= phi <= -30.0 and -77.0 <= psi <= -17.0:
            ss[i] = "H"
        elif -160.0 <= phi <= -50.0 and 90.0 <= psi <= 180.0:
            ss[i] = "E"
        elif -110.0 <= phi <= -40.0 and 20.0 <= psi <= 80.0:
            ss[i] = "T"
        else:
            ss[i] = "C"
    return "".join(ss)


def _pick_files(work: Path):
    tpr = None
    for name in ("md_0_1.tpr", "md.tpr"):
        if (work / name).exists():
            tpr = work / name
            break
    traj = None
    for name in ("md_fit.xtc", "md_0_1.xtc", "md.xtc", "md_noPBC.xtc"):
        if (work / name).exists():
            traj = work / name
            break
    return tpr, traj


def _select_peptide(u, nres_hint: int = 7):
    import MDAnalysis as mda  # noqa: F401
    pep = u.select_atoms("segid B or chainID B")
    if pep.n_residues >= 2:
        return pep
    prot = u.select_atoms("protein")
    rids = list(prot.residues.resids)
    if not rids:
        return pep
    if 531 in rids and 537 in rids:
        return u.select_atoms("resid 531-537")
    last = rids[-nres_hint:]
    return u.select_atoms(f"resid {last[0]}-{last[-1]}")


def compute_from_trajectory(work: Path, nres_hint: int = 7):
    try:
        import MDAnalysis as mda
    except ImportError:
        print(">> [DSSP fallback] MDAnalysis 未安装，无法从轨迹直接计算二级结构。")
        return None, [], None

    tpr, traj = _pick_files(work)
    if tpr is None or traj is None:
        print(">> [DSSP fallback] 未找到 md_0_1.tpr/md.tpr 与轨迹，跳过轨迹计算。")
        return None, [], None

    print(f">> [DSSP fallback] 从轨迹计算 Kabsch-Sander / Ramachandran: {tpr.name} + {traj.name}")
    try:
        u = mda.Universe(str(tpr), str(traj))
    except Exception as exc:
        print(f">> [DSSP fallback] 无法打开轨迹: {exc}")
        return None, [], None

    pep = _select_peptide(u, nres_hint)
    nres = pep.n_residues
    if nres < 3:
        print(">> [DSSP fallback] 肽残基数 < 3，放弃。")
        return None, [], None
    print(f">> [DSSP fallback] 肽残基: {pep.residues.resids[0]}-{pep.residues.resids[-1]} (n={nres})")

    times: List[float] = []
    seqs: List[str] = []
    rama: List[str] = []
    res_ag = list(pep.residues)
    for ts in u.trajectory:
        bb = {name: np.zeros((nres, 3)) for name in ("N", "CA", "C", "O")}
        ok = True
        for i, res in enumerate(res_ag):
            for name in bb:
                ag = res.atoms.select_atoms(f"name {name}")
                if len(ag) == 0:
                    ok = False
                    break
                bb[name][i] = ag.positions[0]
            if not ok:
                break
        if not ok:
            seqs.append("C" * nres)
            rama.append("C" * nres)
        else:
            seqs.append(kabsch_sander_assign(bb))
            rama.append(rama_assign(bb))
        # MDAnalysis frame time is in ps
        times.append(float(ts.time) * 0.001)
    t = np.asarray(times, dtype=float)
    if t.max() <= 0 and len(t) > 1:
        t = infer_times(len(t), 0.1)
    return t, seqs, rama


def try_gmx_dssp(work: Path) -> bool:
    gmx = shutil.which("gmx") or shutil.which("gmx.exe")
    if not gmx:
        return False
    tpr, traj = _pick_files(work)
    if tpr is None or traj is None:
        return False
    ndx = work / "index.ndx"
    sels = ['group "Peptide"', "group Peptide", "Peptide", "resid 531 to 537"]
    for sel in sels:
        cmd = [
            gmx, "dssp",
            "-s", tpr.name, "-f", traj.name,
            "-o", "ss_pep.dat", "-num", "ss_pep_num.xvg",
            "-hmode", "dssp", "-clear",
            "-sel", sel,
        ]
        if ndx.exists():
            cmd[4:4] = ["-n", "index.ndx"]
        print(">> [gmx dssp]", " ".join(cmd))
        try:
            proc = subprocess.run(
                cmd, cwd=str(work), capture_output=True, text=True, timeout=900
            )
        except Exception as exc:
            print(f"   调用失败: {exc}")
            continue
        n_num = 0
        num = work / "ss_pep_num.xvg"
        if num.exists():
            with open(num, encoding="utf-8", errors="ignore") as fh:
                n_num = sum(1 for ln in fh if ln.strip() and not ln.lstrip().startswith(("#", "@")))
        n_dat = 0
        dat = work / "ss_pep.dat"
        if dat.exists():
            with open(dat, encoding="utf-8", errors="ignore") as fh:
                n_dat = sum(1 for ln in fh if ln.strip() and not ln.lstrip().startswith(("#", "@", ";")))
        print(f"   return={proc.returncode}  num_rows={n_num}  dat_rows={n_dat}")
        if proc.returncode == 0 and (n_num >= 10 or n_dat >= 10):
            return True
    return False


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

def write_outputs(
    work: Path,
    times: np.ndarray,
    fracs: Dict[str, np.ndarray],
    seqs: Sequence[str],
    window_ns: float,
    rama_seqs: Optional[Sequence[str]] = None,
    source: str = "unknown",
    prefix: str = "ss_pep",
    title: str = "Peptide",
) -> None:
    perres = work / f"{prefix}_perres.dat"
    with open(perres, "w", encoding="utf-8") as fh:
        fh.write("# time_ns ss_string\n")
        if seqs:
            for t, s in zip(times, seqs):
                fh.write(f"{t:.4f} {s}\n")
    print(f">> [SAVED] {perres}  ({len(seqs)} frames)")

    frac_path = work / f"{prefix}_frac.xvg"
    with open(frac_path, "w", encoding="utf-8") as fh:
        fh.write("# time_ns helix turn bend sheet coil ppii break\n")
        for i, t in enumerate(times):
            fh.write(
                f"{t:.4f} {fracs['helix'][i]:.4f} {fracs['turn'][i]:.4f} "
                f"{fracs['bend'][i]:.4f} {fracs['sheet'][i]:.4f} "
                f"{fracs['coil'][i]:.4f} {fracs['ppii'][i]:.4f} "
                f"{fracs['break'][i]:.4f}\n"
            )
    print(f">> [SAVED] {frac_path}  ({len(times)} frames, per-frame DSSP)")

    bins = bin_fracs(times, fracs, window_ns)
    bins_path = work / f"{prefix}_bins.dat"
    with open(bins_path, "w", encoding="utf-8") as fh:
        fh.write("time_ns helix_frac turn_frac bend_frac coil_frac sheet_frac ppii_frac break_frac n_frames\n")
        for rec in bins:
            fh.write(
                f"{rec['time_ns']:.3f} {rec['helix']:.4f} {rec['turn']:.4f} "
                f"{rec['bend']:.4f} {rec['coil']:.4f} {rec['sheet']:.4f} "
                f"{rec['ppii']:.4f} {rec['break']:.4f} {rec['n_frames']}\n"
            )
    print(f">> [SAVED] {bins_path}  ({len(bins)} windows of {window_ns:g} ns)")

    if rama_seqs:
        rama_fracs = seqs_to_fracs(rama_seqs)
        rama_path = work / f"{prefix}_rama.dat"
        with open(rama_path, "w", encoding="utf-8") as fh:
            fh.write("# time_ns helix_like extended turn_like coil ss_string\n")
            for i, t in enumerate(times):
                fh.write(
                    f"{t:.4f} {rama_fracs['helix'][i]:.4f} {rama_fracs['sheet'][i]:.4f} "
                    f"{rama_fracs['turn'][i]:.4f} {rama_fracs['coil'][i]:.4f} {rama_seqs[i]}\n"
                )
        print(f">> [SAVED] {rama_path}")

    # summary
    n = len(times)
    span = float(times[-1] - times[0]) if n else 0.0
    break_mean = float(np.mean(fracs["break"])) if n else 0.0
    lines = [
        "============================================================",
        f" {title} secondary structure (DSSP) report",
        "============================================================",
        f"source            : {source}",
        f"n_frames          : {n}",
        f"time_range        : {times[0]:.3f} – {times[-1]:.3f} ns" if n else "time_range        : n/a",
        f"bin_window        : {window_ns:g} ns  -> {len(bins)} bins",
        f"SS_string_length  : {Counter(len(s) for s in seqs).most_common(1)[0] if seqs else 'n/a'}",
        "------------------------------------------------------------",
        "full-trajectory fraction (mean ± std):",
    ]
    for k in CATEGORIES:
        if n:
            mu = float(np.mean(fracs[k]) * 100.0)
            sd = float(np.std(fracs[k], ddof=1) * 100.0) if n > 1 else 0.0
            lines.append(f"  {k:<8} {mu:6.2f} ± {sd:5.2f} %")
    if seqs:
        # consensus string
        L = max(len(s) for s in seqs)
        cons = []
        for j in range(L):
            col = [s[j] if j < len(s) else "C" for s in seqs]
            cons.append(Counter(col).most_common(1)[0][0])
        cons_s = "".join(cons)
        if L <= 20:
            lines.append(f"consensus SS      : {cons_s}")
        else:
            n_h = cons_s.count("H") + cons_s.count("G") + cons_s.count("I")
            n_e = cons_s.count("E") + cons_s.count("B")
            lines.append(f"consensus length  : {L} residues  helix_sites={n_h}  sheet_sites={n_e}")
        lines.append("most-frequent DSSP code per residue (first/last 8 if long):")
        show = range(L) if L <= 16 else list(range(8)) + list(range(L - 8, L))
        for j in show:
            col = [s[j] if j < len(s) else "C" for s in seqs]
            code, cnt = Counter(col).most_common(1)[0]
            lines.append(f"  res {j+1:4d}  {code}  ({100.0*cnt/len(col):5.1f} %)")
    if break_mean > 0.05:
        lines.append("------------------------------------------------------------")
        lines.append(
            f"!! 警告: DSSP 将 {break_mean*100:.1f}% 的残基标为链断裂 (=)。"
        )
        lines.append("   4ey6 受体本身有两处缺失环，是三段式多肽；肽分析必须只取 Peptide 组。")
        lines.append("   若断裂出现在 7 肽内部，请确认已使用 md_fit.xtc（去 PBC）再算 DSSP。")
    lines.append("============================================================")
    summary = work / "ss_pep_summary.txt"
    text = "\n".join(lines) + "\n"
    summary.write_text(text, encoding="utf-8")
    print(text)
    print(f">> [SAVED] {summary}")


def _fracs_look_empty(fracs: Dict[str, np.ndarray]) -> bool:
    if not fracs or len(next(iter(fracs.values()))) < 5:
        return True
    stacked = np.vstack([fracs[k] for k in CATEGORIES])
    return float(np.nanmax(stacked)) < 1e-12


def _process_one(work: Path, dat_names, num_name: str, nres: int, trim: bool,
                 prefix: str, title: str, window_ns: float, testing: bool,
                 allow_gmx: bool, gmx_sels: Optional[Sequence[str]],
                 out_dat: str, out_num: str, allow_ks: bool) -> bool:
    times = None
    seqs: List[str] = []
    fracs = None
    rama = None
    source = ""

    for cand in dat_names:
        p = work / cand if not isinstance(cand, Path) else cand
        t, s = parse_ss_dat(p, nres=nres, trim=trim)
        if s and len(s) >= 5:
            times, seqs = t, s
            fracs = seqs_to_fracs(seqs)
            source = f"parsed {p.name} ({len(seqs)} frames, L={len(seqs[0])})"
            print(f">> {title}: {p.name}  frames={len(seqs)}  SS_len={len(seqs[0])}")
            break

    t_num, f_num = parse_ss_num_xvg(work / num_name, nres=max(nres, 50) if not trim else nres)
    if f_num is not None and (fracs is None or _fracs_look_empty(fracs)):
        times, fracs = t_num, f_num
        source = f"parsed {num_name} ({len(times)} frames)"
        print(f">> {title}: {num_name}  frames={len(times)}")

    need = times is None or fracs is None or _fracs_look_empty(fracs) or (len(times) < 10)
    # complex must be long
    if not trim and seqs and len(seqs[0]) < 50:
        print(f">> {title}: SS length {len(seqs[0])} is too short for a complex; will retry gmx")
        need = True

    if need and allow_gmx and gmx_sels:
        print(f">> {title}: running gmx dssp ...")
        if try_gmx_dssp(work, sels=gmx_sels, out_dat=out_dat, out_num=out_num):
            t, s = parse_ss_dat(work / out_dat, nres=nres, trim=trim)
            t2, f2 = parse_ss_num_xvg(work / out_num, nres=max(nres, 50) if not trim else nres)
            if s and len(s) >= 5 and (trim or len(s[0]) >= 50):
                times, seqs, fracs = t, s, seqs_to_fracs(s)
                source = f"gmx dssp -o {out_dat}"
                need = False
            elif f2 is not None:
                times, fracs = t2, f2
                source = f"gmx dssp -num {out_num}"
                need = False

    if need and allow_ks and trim:
        t3, s3, r3 = compute_from_trajectory(work, nres_hint=nres)
        if s3 and len(s3) >= 5:
            times, seqs, rama = t3, s3, r3
            fracs = seqs_to_fracs(seqs)
            source = "Kabsch-Sander fallback"
            need = False

    if need or times is None or fracs is None:
        print(f"!!! {title} DSSP failed")
        return False

    if window_ns <= 0:
        window_ns = auto_window_ns(float(times[0]), float(times[-1]), testing=testing)
    write_outputs(work, times, fracs, seqs, window_ns, rama_seqs=rama,
                  source=source, prefix=prefix, title=title)
    return True


def run(work: Path, nres: int, window_ns: float, testing: bool, allow_gmx: bool,
        target: str = "both") -> int:
    work = work.resolve()
    print("=" * 60)
    print(f" DSSP rebuild (target={target}) in: {work}")
    print("=" * 60)

    ok_c = ok_p = True
    if target in ("complex", "both"):
        ok_c = _process_one(
            work,
            dat_names=["ss_complex.dat", "ss_ache.dat"],
            num_name="ss_complex_num.xvg",
            nres=537, trim=False,
            prefix="ss_complex", title="Complex(AChE+Peptide)",
            window_ns=window_ns, testing=testing, allow_gmx=allow_gmx,
            gmx_sels=[
                'group "Protein" or group "Peptide"',
                "group Protein or group Peptide",
                "group Protein",
                'group "AChE" or group "Peptide"',
                "resid 1 to 537",
            ],
            out_dat="ss_complex.dat", out_num="ss_complex_num.xvg",
            allow_ks=False,
        )
    if target in ("peptide", "both"):
        ok_p = _process_one(
            work,
            dat_names=["ss_pep.dat", "ss_pep.sc"],
            num_name="ss_pep_num.xvg",
            nres=nres, trim=True,
            prefix="ss_pep", title="Peptide",
            window_ns=window_ns, testing=testing, allow_gmx=allow_gmx,
            gmx_sels=['group "Peptide"', "group Peptide", "Peptide", "resid 531 to 537"],
            out_dat="ss_pep.dat", out_num="ss_pep_num.xvg",
            allow_ks=True,
        )
    if target == "complex":
        return 0 if ok_c else 1
    if target == "peptide":
        return 0 if ok_p else 1
    # both: succeed if complex worked (main request); peptide is extra
    return 0 if ok_c or ok_p else 1


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Rebuild per-frame DSSP tables")
    ap.add_argument("-d", "--dir", default=".", help="MD working directory")
    ap.add_argument("--nres", type=int, default=7, help="Peptide residue count")
    ap.add_argument("--window-ns", type=float, default=0.0, help="Bin width; 0 = auto")
    ap.add_argument("--target", choices=("complex", "peptide", "both"), default="both")
    ap.add_argument("--testing", action="store_true")
    ap.add_argument("--no-gmx", action="store_true", help="Do not invoke gmx dssp")
    args = ap.parse_args(argv)
    testing = args.testing or (str(__import__("os").environ.get("TESTING", "0")) == "1")
    return run(Path(args.dir), args.nres, args.window_ns, testing,
               allow_gmx=not args.no_gmx, target=args.target)


if __name__ == "__main__":
    sys.exit(main())
