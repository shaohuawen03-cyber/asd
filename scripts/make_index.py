#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_index.py
=============

Build an index file with the groups the analysis stage needs:

    Receptor            AChE  (all protein moleculetypes except the peptide)
    Peptide             the docked peptide
    Complex             Receptor + Peptide
    Receptor_backbone   N, CA, C of the receptor
    Peptide_backbone    N, CA, C of the peptide
    Complex_backbone    N, CA, C of both
    Water_and_ions      SOL + NA + CL

Chain IDs do not survive into the .gro file, so the peptide is located the
reliable way: the ``[ molecules ]`` order in ``topol.top`` matches the atom
order in the coordinate file, so summing the atom counts of the preceding
moleculetypes gives the peptide's atom range exactly.

The default GROMACS groups produced by ``gmx make_ndx`` are kept and the custom
groups are appended, so anything a normal workflow expects is still there.

Usage
-----
    python3 scripts/make_index.py --struct solv_ions.gro --top topol.top \\
            --prepared-json prepared.json -o index.ndx
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

BACKBONE_NAMES = {"N", "CA", "C"}


# --------------------------------------------------------------------------- #
#  Topology parsing
# --------------------------------------------------------------------------- #

def strip_comment(line: str) -> str:
    return line.split(";", 1)[0].rstrip()


def read_top_with_includes(path: str, _depth: int = 0) -> List[Tuple[str, str]]:
    """
    Flatten a .top and its local #includes into [(source_file, line), ...].
    Force-field includes (paths containing '.ff/') are skipped - we only care
    about moleculetype definitions written by pdb2gmx.
    """
    if _depth > 10:
        return []
    out: List[Tuple[str, str]] = []
    base = os.path.dirname(os.path.abspath(path))
    with open(path, "r", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            m = re.match(r'\s*#include\s+"([^"]+)"', line)
            if m:
                inc = m.group(1)
                if ".ff/" in inc.replace("\\", "/"):
                    continue
                inc_path = inc if os.path.isabs(inc) else os.path.join(base, inc)
                if os.path.isfile(inc_path):
                    out.extend(read_top_with_includes(inc_path, _depth + 1))
                continue
            out.append((path, line))
    return out


def parse_topology(top_path: str):
    """
    Return (moltype_atoms, molecule_list) where

      moltype_atoms : {moleculetype_name: [(atom_name, residue_name), ...]}
      molecule_list : [(moleculetype_name, count), ...]  in [ molecules ] order
    """
    lines = read_top_with_includes(top_path)

    moltype_atoms: Dict[str, List[Tuple[str, str]]] = {}
    molecules: List[Tuple[str, int]] = []

    section = None
    current_mol: Optional[str] = None

    for _src, raw in lines:
        line = strip_comment(raw).strip()
        if not line:
            continue
        if line.startswith("#"):        # #ifdef / #endif / #define ...
            continue
        m = re.match(r"\[\s*(.+?)\s*\]", line)
        if m:
            section = m.group(1).lower()
            if section == "moleculetype":
                current_mol = None
            continue

        if section == "moleculetype":
            current_mol = line.split()[0]
            moltype_atoms.setdefault(current_mol, [])
        elif section == "atoms" and current_mol:
            f = line.split()
            # nr type resnr residue atom cgnr charge [mass]
            if len(f) >= 5:
                moltype_atoms[current_mol].append((f[4], f[3]))
        elif section == "molecules":
            f = line.split()
            if len(f) >= 2:
                try:
                    molecules.append((f[0], int(f[1])))
                except ValueError:
                    pass

    return moltype_atoms, molecules


# --------------------------------------------------------------------------- #
#  Index writing
# --------------------------------------------------------------------------- #

def run_default_make_ndx(gmx: str, struct: str, out: str) -> bool:
    """
    Ask make_ndx for the standard groups.  A bare 'q' makes some GROMACS
    versions exit without writing anything, so an explicit no-op selection
    ('0' = System) is issued first to mark the set as modified.
    """
    try:
        subprocess.run([gmx, "make_ndx", "-f", struct, "-o", out],
                       input="0\nq\n", text=True, check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"WARNING: default make_ndx failed ({exc}); "
              f"writing custom groups only", file=sys.stderr)
        return False
    if not os.path.isfile(out):
        print("WARNING: make_ndx produced no file; writing custom groups only",
              file=sys.stderr)
        return False
    return True


def read_ndx(path: str) -> "list[tuple[str, list[int]]]":
    groups: List[Tuple[str, List[int]]] = []
    name = None
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
                idx.extend(int(t) for t in line.split() if t.lstrip("-").isdigit())
    if name is not None:
        groups.append((name, idx))
    return groups


def write_ndx(path: str, groups: "list[tuple[str, list[int]]]") -> None:
    with open(path, "w") as fh:
        for name, idx in groups:
            fh.write(f"[ {name} ]\n")
            for i in range(0, len(idx), 15):
                fh.write(" ".join(f"{v:6d}" for v in idx[i:i + 15]) + "\n")
            if not idx:
                fh.write("\n")


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Index groups for AChE-peptide MD")
    ap.add_argument("--struct", required=True, help="solv_ions.gro / .tpr")
    ap.add_argument("--top", required=True, help="topol.top")
    ap.add_argument("--prepared-json", default=None,
                    help="report from prepare_structure.py (for peptide chain)")
    ap.add_argument("--peptide-moltypes", default=None,
                    help="comma-separated moleculetype names forming the peptide "
                         "(overrides autodetection)")
    ap.add_argument("--gmx", default="gmx")
    ap.add_argument("-o", "--output", default="index.ndx")
    args = ap.parse_args(argv)

    moltype_atoms, molecules = parse_topology(args.top)
    if not molecules:
        print("ERROR: no [ molecules ] section found in topology", file=sys.stderr)
        return 2

    # ---- which moleculetypes are the peptide? ------------------------------
    protein_mols = [(n, c) for n, c in molecules
                    if n in moltype_atoms and moltype_atoms[n]]
    if not protein_mols:
        print("ERROR: no protein moleculetypes found", file=sys.stderr)
        return 2

    if args.peptide_moltypes:
        peptide_names = {s.strip() for s in args.peptide_moltypes.split(",")}
    else:
        peptide_names = None
        if args.prepared_json and os.path.isfile(args.prepared_json):
            rep = json.load(open(args.prepared_json))
            pep_chain = rep.get("peptide_chain")
            segs = rep.get("segments", [])
            if pep_chain:
                # position of the peptide segments among all segments; pdb2gmx
                # emits moleculetypes in exactly that order
                pep_positions = [i for i, s in enumerate(segs)
                                 if s.get("chain") == pep_chain]
                if pep_positions and len(segs) == len(protein_mols):
                    peptide_names = {protein_mols[i][0] for i in pep_positions}
        if peptide_names is None:
            # fallback: the smallest protein moleculetype is the peptide
            smallest = min(protein_mols, key=lambda t: len(moltype_atoms[t[0]]))
            peptide_names = {smallest[0]}
            print(f"NOTE: peptide autodetected as smallest moleculetype "
                  f"'{smallest[0]}' ({len(moltype_atoms[smallest[0]])} atoms)")

    # ---- walk the coordinate file atom order -------------------------------
    receptor: List[int] = []
    peptide: List[int] = []
    receptor_bb: List[int] = []
    peptide_bb: List[int] = []
    water_ions: List[int] = []

    serial = 0
    for name, count in molecules:
        atoms = moltype_atoms.get(name, [])
        if not atoms:
            # solvent / ions: definition lives in the force-field include
            n_per = {"SOL": 3, "NA": 1, "CL": 1, "K": 1, "MG": 1, "CA": 1,
                     "ZN": 1, "TIP3": 3, "HOH": 3, "WAT": 3}.get(name.upper())
            if n_per is None:
                print(f"ERROR: unknown molecule '{name}' with no [ atoms ] "
                      f"section; pass --peptide-moltypes and check the topology",
                      file=sys.stderr)
                return 2
            total = n_per * count
            water_ions.extend(range(serial + 1, serial + total + 1))
            serial += total
            continue

        is_pep = name in peptide_names
        for _ in range(count):
            for atom_name, _resname in atoms:
                serial += 1
                if is_pep:
                    peptide.append(serial)
                    if atom_name in BACKBONE_NAMES:
                        peptide_bb.append(serial)
                else:
                    receptor.append(serial)
                    if atom_name in BACKBONE_NAMES:
                        receptor_bb.append(serial)

    complex_all = sorted(receptor + peptide)
    complex_bb = sorted(receptor_bb + peptide_bb)

    custom = [
        ("Receptor", receptor),
        ("Peptide", peptide),
        ("Complex", complex_all),
        ("Receptor_backbone", receptor_bb),
        ("Peptide_backbone", peptide_bb),
        ("Complex_backbone", complex_bb),
    ]
    if water_ions:
        custom.append(("Water_and_ions", water_ions))

    # ---- merge with the default groups -------------------------------------
    # NB: GROMACS insists on a .ndx extension for the output of make_ndx
    tmp = os.path.splitext(args.output)[0] + "_default.ndx"
    groups: List[Tuple[str, List[int]]] = []
    if run_default_make_ndx(args.gmx, args.struct, tmp):
        existing = read_ndx(tmp)
        # the no-op '0' selection duplicates "System"; keep the first of each
        deduped: List[Tuple[str, List[int]]] = []
        taken = set()
        for n, i in existing:
            if n in taken:
                continue
            taken.add(n)
            deduped.append((n, i))
        groups = deduped + [(n, i) for n, i in custom if n not in taken]
        os.remove(tmp)
    else:
        groups = custom

    write_ndx(args.output, groups)

    print("-" * 60)
    print(f"index written: {args.output}")
    for n, i in custom:
        print(f"   {n:<20s} {len(i):8d} atoms")
    print(f"   {'TOTAL parsed':<20s} {serial:8d} atoms")
    print("-" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
