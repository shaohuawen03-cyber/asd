#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
water_bridges.py
================

Reproduces Figure 6 / Table 2 of the reference paper:

    "A bridging water molecule is a water molecule that is bound to both Abeta
     and AChE in at least one frame of the 5000 frames of the production
     dynamics.  We assumed that each bridge is a single continuous interaction
     and that a given water molecule can form several bridging interactions.
     The distance cut-off for the hydrogen bonds was 3 A and the angle cut-off
     was 135 degrees."

For every peptide residue the script reports

  * ``n_bridge_waters``   distinct water molecules that ever bridged it
  * ``n_bridges``         number of *continuous* bridging interactions
                          (a water leaving and coming back counts twice - this
                          is why the paper's bridge count exceeds the molecule
                          count)
  * lifetime statistics of those interactions, including the fraction that are
    short-lived (a single frame)

Definition used
---------------
A water W bridges peptide residue R and the receptor in a frame when, in that
same frame,

    min distance(W, atoms of R)        <= cutoff   AND
    min distance(W, receptor atoms)    <= cutoff

with ``cutoff`` = 0.35 nm by default. The paper's 3 A criterion is a
donor-acceptor *hydrogen-bond* distance; here the geometric criterion is
applied on a heavy-atom basis, which is the standard GROMACS convention
(``--cutoff 0.3`` reproduces the literal 3 A value).

Implementation
--------------
Water selections per frame come from ``gmx select -on``, which writes one index
group per frame; those groups are intersected in Python.  Only the standard
library is needed.

Usage
-----
    python3 scripts/water_bridges.py -s md.tpr -f md_fit.xtc -n index.ndx \\
            -o analysis/ --stride 5
"""

from __future__ import annotations

import argparse
import collections
import os
import re
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional, Sequence, Tuple


# --------------------------------------------------------------------------- #
#  helpers
# --------------------------------------------------------------------------- #

def run_gmx(cmd: List[str], check: bool = True) -> str:
    proc = subprocess.run(cmd, input="", text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if check and proc.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(cmd)}\n{proc.stdout}")
    return proc.stdout


def read_ndx_ordered(path: str) -> List[Tuple[str, List[int]]]:
    """Read an ndx preserving order and duplicate names (one group per frame)."""
    groups: List[Tuple[str, List[int]]] = []
    name: Optional[str] = None
    idx: List[int] = []
    with open(path) as fh:
        for line in fh:
            m = re.match(r"\s*\[\s*(.+?)\s*\]", line)
            if m:
                if name is not None:
                    groups.append((name, idx))
                name = m.group(1)
                idx = []
            elif name is not None:
                idx.extend(int(t) for t in line.split()
                           if t.lstrip("-").isdigit())
    if name is not None:
        groups.append((name, idx))
    return groups


def read_ndx_dict(path: str) -> Dict[str, List[int]]:
    out: Dict[str, List[int]] = {}
    for n, i in read_ndx_ordered(path):
        out.setdefault(n, []).extend(i)
    return out


def atom_table(gmx: str, tpr: str) -> Dict[int, Tuple[int, str, str]]:
    """{global_index_1based: (resnum, resname, atomname)} in coordinate order."""
    table: Dict[int, Tuple[int, str, str]] = {}
    with tempfile.TemporaryDirectory() as td:
        pdb = os.path.join(td, "topology.pdb")
        run_gmx([gmx, "editconf", "-f", tpr, "-o", pdb])
        serial = 0
        with open(pdb, errors="replace") as fh:
            for line in fh:
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                serial += 1
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    resnum = 0
                table[serial] = (resnum, line[17:20].strip(),
                                 line[12:16].strip())
    return table


def traj_timing(gmx: str, traj: str) -> Tuple[int, float]:
    """
    (number_of_frames, frame_spacing_ps) from `gmx check`, whose summary looks
    like::

        Item        #frames Timestep (ps)
        Step            11    0.02

    The timestep is a float, so it must not be matched as an integer.
    """
    out = run_gmx([gmx, "check", "-f", traj], check=False)
    n_frames, dt = 0, 0.0
    m = re.search(r"^\s*Step\s+(\d+)\s+([0-9.eE+-]+)", out, re.MULTILINE)
    if m:
        n_frames = int(m.group(1))
        try:
            dt = float(m.group(2))
        except ValueError:
            dt = 0.0
    if dt == 0.0:
        # single-frame trajectory, or an unparsable header: fall back to the
        # first two "Reading frame" lines
        times = [float(t) for t in
                 re.findall(r"Reading frame\s+\d+\s+time\s+([0-9.eE+-]+)", out)]
        if len(times) >= 2:
            dt = times[1] - times[0]
        n_frames = n_frames or len(times)
    return n_frames, dt


# --------------------------------------------------------------------------- #
#  Core analysis
# --------------------------------------------------------------------------- #

def bridging_waters_per_frame(gmx: str, tpr: str, traj: str, ndx: str,
                              sel: str, stride: int,
                              begin: Optional[float], end: Optional[float],
                              workdir: str, tag: str) -> List[List[int]]:
    """Run one gmx select and return the selected atom indices for each frame."""
    out_ndx = os.path.join(workdir, f"sel_{tag}.ndx")
    # -selrpos atom: distances are measured from the water oxygen itself, not
    # from a residue centre of mass, which is what a hydrogen-bond style
    # criterion requires.
    cmd = [gmx, "select", "-s", tpr, "-f", traj, "-n", ndx,
           "-select", sel, "-on", out_ndx, "-selrpos", "atom"]
    if stride and stride > 0:
        cmd += ["-dt", str(stride)]
    if begin is not None:
        cmd += ["-b", str(begin)]
    if end is not None:
        cmd += ["-e", str(end)]
    run_gmx(cmd)
    return [idx for _n, idx in read_ndx_ordered(out_ndx)]


def water_id(atom_index: int, atoms: Dict[int, Tuple[int, str, str]]) -> int:
    """Residue number of the water an atom belongs to (its molecule id)."""
    return atoms.get(atom_index, (0, "", ""))[0]


def summarise_lifetimes(presence: Dict[int, List[int]], n_frames: int,
                        dt_ns: float):
    """
    presence: {water_resid: [frame indices where it bridges]}
    Returns (n_waters, n_bridges, lifetimes_ns list).
    """
    lifetimes: List[float] = []
    n_bridges = 0
    for _wid, frames in presence.items():
        frames = sorted(frames)
        run_len = 1
        for a, b in zip(frames, frames[1:]):
            if b == a + 1:
                run_len += 1
            else:
                n_bridges += 1
                lifetimes.append(run_len * dt_ns)
                run_len = 1
        n_bridges += 1
        lifetimes.append(run_len * dt_ns)
    return len(presence), n_bridges, lifetimes


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Water-mediated bridges between peptide and receptor")
    ap.add_argument("-s", "--tpr", required=True)
    ap.add_argument("-f", "--traj", required=True)
    ap.add_argument("-n", "--ndx", required=True)
    ap.add_argument("-o", "--outdir", default="analysis")
    ap.add_argument("--cutoff", type=float, default=0.35,
                    help="heavy-atom distance cut-off in nm "
                         "(0.35 default; use 0.3 for the paper's literal 3 A)")
    ap.add_argument("--dt", dest="stride", type=float, default=0.0,
                    help="only analyse frames every DT picoseconds "
                         "(0 = use every frame; passed to gmx as -dt)")
    ap.add_argument("--peptide-group", default="Peptide")
    ap.add_argument("--receptor-group", default="Receptor")
    ap.add_argument("--water-select", default="resname SOL and name OW",
                    help="selection identifying water oxygens")
    ap.add_argument("-b", "--begin", type=float, default=None)
    ap.add_argument("-e", "--end", type=float, default=None)
    ap.add_argument("--gmx", default="gmx")
    args = ap.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    print("=" * 66)
    print("water_bridges.py  (paper Fig 6 / Table 2)")
    print(f"  cutoff {args.cutoff} nm ({args.cutoff*10:.1f} A), "
          f"dt {args.stride or 'all frames'}")
    print("=" * 66)

    groups = read_ndx_dict(args.ndx)
    for g in (args.peptide_group, args.receptor_group):
        if g not in groups:
            print(f"ERROR: group '{g}' not in {args.ndx}", file=sys.stderr)
            return 2

    print("  reading topology ...")
    atoms = atom_table(args.gmx, args.tpr)

    pep_atoms = groups[args.peptide_group]
    pep_residues: "collections.OrderedDict[Tuple[int,str], List[int]]" = \
        collections.OrderedDict()
    for a in pep_atoms:
        if a in atoms:
            rnum, rname, _ = atoms[a]
            pep_residues.setdefault((rnum, rname), []).append(a)
    print(f"  peptide: {len(pep_residues)} residues")

    cut = args.cutoff
    rec = args.receptor_group
    wsel = args.water_select

    results = []
    with tempfile.TemporaryDirectory() as td:
        # waters close to the receptor, frame by frame (computed once)
        print("  locating waters near the receptor ...")
        near_rec = bridging_waters_per_frame(
            args.gmx, args.tpr, args.traj, args.ndx,
            f'({wsel}) and within {cut} of group "{rec}"',
            args.stride, args.begin, args.end, td, "rec")
        n_frames = len(near_rec)
        if n_frames == 0:
            print("ERROR: no frames selected", file=sys.stderr)
            return 2
        near_rec_sets = [set(water_id(a, atoms) for a in fr) for fr in near_rec]
        print(f"  {n_frames} frames analysed")

        # spacing between the frames we actually analysed.  --stride is
        # passed to gmx as -dt (picoseconds), so the analysed spacing is the
        # larger of the trajectory's own spacing and the requested one.
        _n_all, dt_traj_ps = traj_timing(args.gmx, args.traj)
        dt_ps = max(dt_traj_ps, float(args.stride)) if dt_traj_ps > 0 \
            else float(args.stride)
        dt_ns = dt_ps / 1000.0
        print(f"  frame spacing analysed: {dt_ps:g} ps ({dt_ns:g} ns)")

        for i, ((rnum, rname), aidx) in enumerate(pep_residues.items(), 1):
            tag = f"{rname}{rnum}"
            idx_expr = " ".join(str(a) for a in aidx)
            per_frame = bridging_waters_per_frame(
                args.gmx, args.tpr, args.traj, args.ndx,
                f'({wsel}) and within {cut} of (atomnr {idx_expr})',
                args.stride, args.begin, args.end, td, f"r{rnum}")

            presence: Dict[int, List[int]] = collections.defaultdict(list)
            for fi, fr in enumerate(per_frame[:n_frames]):
                near_r = set(water_id(a, atoms) for a in fr)
                for wid in near_r & near_rec_sets[fi]:
                    presence[wid].append(fi)

            n_w, n_b, lifetimes = summarise_lifetimes(presence, n_frames, dt_ns)
            short = sum(1 for l in lifetimes if l <= dt_ns + 1e-9)
            results.append({
                "residue": tag,
                "n_bridge_waters": n_w,
                "n_bridges": n_b,
                "max_lifetime_ns": max(lifetimes) if lifetimes else 0.0,
                "mean_lifetime_ns": (sum(lifetimes) / len(lifetimes)
                                     if lifetimes else 0.0),
                "frac_single_frame": short / n_b if n_b else 0.0,
            })
            print(f"    [{i:3d}/{len(pep_residues)}] {tag:<8s} "
                  f"waters {n_w:5d}  bridges {n_b:5d}  "
                  f"max life {results[-1]['max_lifetime_ns']*1000:9.3f} ps")

    out = os.path.join(args.outdir, "water_bridges.dat")
    with open(out, "w") as fh:
        fh.write("# Water-mediated bridges between the peptide and the "
                 "receptor (paper Fig 6 / Table 2)\n")
        fh.write(f"# cutoff {args.cutoff} nm, {n_frames} frames analysed, "
                 f"frame spacing {dt_ns:.4f} ns\n")
        fh.write("# residue  n_bridge_waters  n_bridges  max_lifetime_ns  "
                 "mean_lifetime_ns  frac_single_frame\n")
        fh.write("# 'frac_single_frame' = fraction of bridges lasting one "
                 "frame only (paper: most bridges are short-lived)\n")
        for r in results:
            fh.write(f"{r['residue']:<10s} {r['n_bridge_waters']:10d} "
                     f"{r['n_bridges']:12d} {r['max_lifetime_ns']:14.4f} "
                     f"{r['mean_lifetime_ns']:16.4f} "
                     f"{r['frac_single_frame']:16.3f}\n")

    ranked = sorted(results, key=lambda r: r["n_bridges"], reverse=True)
    print()
    print("  " + "-" * 60)
    print("  Peptide residues ranked by number of bridging interactions")
    print("  " + "-" * 60)
    for r in ranked[:12]:
        print(f"    {r['residue']:<8s} {r['n_bridges']:6d} bridges, "
              f"{r['n_bridge_waters']:5d} distinct waters, "
              f"longest {r['max_lifetime_ns']*1000:.3f} ps "
              f"({r['frac_single_frame']*100:.0f}% single-frame)")
    print("  " + "-" * 60)
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
