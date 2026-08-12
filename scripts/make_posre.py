#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_posre.py
=============

The paper restrains **heavy atoms of the protein and the peptide** with a
harmonic force constant of **3 kcal/mol/A^2**.

In GROMACS units that is

    3 kcal/mol/A^2 x 4.184 kJ/kcal x 100 A^2/nm^2 = 1255.2 kJ/mol/nm^2

``gmx pdb2gmx`` already writes one ``posre_*.itp`` per moleculetype, but it
hard-codes 1000 kJ/mol/nm^2 and guards them with ``#ifdef POSRES``.  This script

  1. rewrites every ``posre_*.itp`` with the correct 1255.2 kJ/mol/nm^2 and
  2. switches the guard in the topology from ``POSRES`` to ``POSRES_HEAVY``,
     which is the define used by ``mdp/em.mdp``, ``mdp/nvt.mdp`` and
     ``mdp/npt.mdp``.

pdb2gmx restrains exactly the heavy atoms that were present in the input
coordinates (hydrogens it added itself are never restrained), which is precisely
the paper's definition, so the atom selection is reused as-is.

Idempotent: running it twice is harmless.

Usage
-----
    python3 scripts/make_posre.py --dir work/alllhrc
    python3 scripts/make_posre.py --dir work/x --fc-kcal 3.0
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys
from typing import List, Optional

KCAL_TO_KJ = 4.184
A2_PER_NM2 = 100.0


def convert_fc(fc_kcal: float) -> float:
    """kcal/mol/A^2 -> kJ/mol/nm^2"""
    return fc_kcal * KCAL_TO_KJ * A2_PER_NM2


def rewrite_posre(path: str, fc: float, fc_kcal: float) -> int:
    """Rewrite the force constants in one posre itp; returns #atoms restrained."""
    with open(path, "r") as fh:
        lines = fh.readlines()

    out: List[str] = []
    in_section = False
    n = 0
    fc_str = f"{fc:.1f}"

    for raw in lines:
        stripped = raw.strip()
        m = re.match(r"\[\s*(.+?)\s*\]", stripped)
        if m:
            in_section = m.group(1).lower() == "position_restraints"
            out.append(raw)
            continue
        if not in_section or not stripped or stripped.startswith(";"):
            out.append(raw)
            continue
        f = stripped.split()
        if len(f) >= 5:
            out.append(f"{int(f[0]):6d} {int(f[1]):5d} "
                       f"{fc_str:>10s} {fc_str:>10s} {fc_str:>10s}\n")
            n += 1
        else:
            out.append(raw)

    header = (f"; Position restraints on heavy atoms\n"
              f"; force constant = {fc_kcal} kcal/mol/A^2 "
              f"= {fc_str} kJ/mol/nm^2  (reference paper)\n"
              f"; enabled with -DPOSRES_HEAVY\n")
    body = [l for l in out if not l.startswith("; Position restraints on heavy")
            and not l.startswith("; force constant =")
            and not l.startswith("; enabled with -DPOSRES_HEAVY")]

    with open(path, "w") as fh:
        fh.write(header)
        fh.writelines(body)
    return n


def retag_topology(path: str, old: str = "POSRES",
                   new: str = "POSRES_HEAVY") -> int:
    """Change '#ifdef POSRES' -> '#ifdef POSRES_HEAVY' (leaves POSRES_WATER)."""
    with open(path, "r") as fh:
        text = fh.read()
    # word-boundary match so POSRES_WATER / POSRES_HEAVY are untouched
    new_text, count = re.subn(rf"#ifdef\s+{old}\b(?!_)", f"#ifdef {new}", text)
    if count:
        with open(path, "w") as fh:
            fh.write(new_text)
    return count


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Set position restraints to the paper's 3 kcal/mol/A^2")
    ap.add_argument("--dir", required=True, help="system working directory")
    ap.add_argument("--fc-kcal", type=float, default=3.0,
                    help="force constant in kcal/mol/A^2 (default: 3.0)")
    ap.add_argument("--define", default="POSRES_HEAVY",
                    help="preprocessor define to guard the restraints")
    args = ap.parse_args(argv)

    if not os.path.isdir(args.dir):
        print(f"ERROR: no such directory: {args.dir}", file=sys.stderr)
        return 2

    fc = convert_fc(args.fc_kcal)
    posres = sorted(glob.glob(os.path.join(args.dir, "posre*.itp")))
    if not posres:
        print(f"ERROR: no posre*.itp in {args.dir} - run pdb2gmx first",
              file=sys.stderr)
        return 2

    print("=" * 64)
    print("make_posre.py")
    print(f"  force constant : {args.fc_kcal} kcal/mol/A^2 = {fc:.1f} kJ/mol/nm^2")
    print(f"  guard define   : -D{args.define}")
    print("=" * 64)

    total = 0
    for p in posres:
        n = rewrite_posre(p, fc, args.fc_kcal)
        total += n
        print(f"  {os.path.basename(p):<40s} {n:7d} atoms restrained")

    tops = sorted(glob.glob(os.path.join(args.dir, "topol*.top"))
                  + glob.glob(os.path.join(args.dir, "topol*.itp")))
    retagged = 0
    for t in tops:
        c = retag_topology(t, "POSRES", args.define)
        retagged += c
        if c:
            print(f"  {os.path.basename(t):<40s} {c:7d} #ifdef POSRES "
                  f"-> #ifdef {args.define}")

    print("-" * 64)
    print(f"  total restrained atoms : {total}")
    if retagged == 0:
        print(f"  (no #ifdef POSRES found - already tagged {args.define}?)")
    print("=" * 64)
    return 0


if __name__ == "__main__":
    sys.exit(main())
