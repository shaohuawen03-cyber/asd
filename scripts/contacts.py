#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
contacts.py
===========

Reproduces Figure 5 of the reference paper:

    "All other contacts between a pair of atoms within a distance of 7 A were
     defined as non-native contacts.  The number of contacts of the Abeta
     residues with AChE averaged over 5000 frames is presented in Fig 5A"
    "...the number of intramolecular contacts formed within Abeta (Fig 5B)"

For every peptide residue it reports

  * ``n_contacts``            mean number of receptor atoms within the cut-off
  * ``n_frames_in_contact``   how many frames the residue touches the receptor
  * ``occupancy``             that count as a fraction of all frames
  * the receptor residues it touches most often (the paper's Table 1)

and the same for peptide-internal (intramolecular) contacts.

Distances are computed by ``gmx select`` (which understands PBC and the tpr
topology) and the per-residue book-keeping is done here, so no MDAnalysis /
numpy dependency is needed.

Usage
-----
    python3 scripts/contacts.py -s md.tpr -f md_fit.xtc -n index.ndx \\
            --cutoff 0.7 -o analysis/
"""

from __future__ import annotations

import argparse
import collections
import os
import re
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional, Tuple


# --------------------------------------------------------------------------- #
#  gmx helpers
# --------------------------------------------------------------------------- #

def run_gmx(cmd: List[str], stdin: str = "", check: bool = True) -> str:
    proc = subprocess.run(cmd, input=stdin, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if check and proc.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(cmd)}\n{proc.stdout}")
    return proc.stdout


def parse_xvg(path: str) -> List[List[float]]:
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "@")):
                continue
            f = line.split()
            if f:
                try:
                    rows.append([float(v) for v in f])
                except ValueError:
                    pass
    return rows


# --------------------------------------------------------------------------- #
#  Structure info straight out of the tpr via gmx dump
# --------------------------------------------------------------------------- #

def atom_table(gmx: str, tpr: str) -> Dict[int, Tuple[int, str, str, str]]:
    """
    {global_atom_index_1based: (residue_number, residue_name, atom_name, chain)}

    ``gmx dump`` indexes atoms *per moleculetype*, which does not map onto the
    global numbering used by index files, so the topology is instead dumped to
    a PDB with ``gmx editconf``: that file is written in exactly the global
    atom order and carries residue numbers, names and chain IDs.
    """
    table: Dict[int, Tuple[int, str, str, str]] = {}
    with tempfile.TemporaryDirectory() as td:
        pdb = os.path.join(td, "topology.pdb")
        run_gmx([gmx, "editconf", "-f", tpr, "-o", pdb], check=True)

        serial = 0
        with open(pdb, errors="replace") as fh:
            for line in fh:
                if not line.startswith(("ATOM  ", "HETATM")):
                    continue
                serial += 1
                atom_name = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21]
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    resnum = 0
                table[serial] = (resnum, resname, atom_name, chain)
    return table


def read_ndx_groups(path: str) -> "collections.OrderedDict[str, List[int]]":
    groups: "collections.OrderedDict[str, List[int]]" = collections.OrderedDict()
    name = None
    with open(path) as fh:
        for line in fh:
            m = re.match(r"\s*\[\s*(.+?)\s*\]", line)
            if m:
                name = m.group(1)
                groups.setdefault(name, [])
            elif name is not None:
                groups[name].extend(int(t) for t in line.split()
                                    if t.lstrip("-").isdigit())
    return groups


# --------------------------------------------------------------------------- #
#  Main analysis
# --------------------------------------------------------------------------- #

def analyse(gmx: str, tpr: str, xtc: str, ndx: str, cutoff: float,
            outdir: str, pep_group: str, rec_group: str,
            begin: Optional[float], end: Optional[float]) -> int:

    os.makedirs(outdir, exist_ok=True)
    groups = read_ndx_groups(ndx)
    for g in (pep_group, rec_group):
        if g not in groups:
            print(f"ERROR: group '{g}' not in {ndx}", file=sys.stderr)
            return 2

    print("  reading topology from the tpr ...")
    atoms = atom_table(gmx, tpr)
    if not atoms:
        print("ERROR: could not parse atoms from the tpr", file=sys.stderr)
        return 2

    pep_atoms = groups[pep_group]
    pep_residues: "collections.OrderedDict[Tuple[int,str], List[int]]" = \
        collections.OrderedDict()
    for a in pep_atoms:
        if a not in atoms:
            continue
        rnum, rname, _, chain = atoms[a]
        pep_residues.setdefault((rnum, rname, chain), []).append(a)

    print(f"  peptide: {len(pep_residues)} residues, {len(pep_atoms)} atoms")

    time_args: List[str] = []
    if begin is not None:
        time_args += ["-b", str(begin)]
    if end is not None:
        time_args += ["-e", str(end)]

    # ------------------------------------------------------------------ #
    #  per-residue contact counts, one gmx select call per peptide residue
    # ------------------------------------------------------------------ #
    inter_rows: List[Dict] = []
    intra_rows: List[Dict] = []
    partner_counts: Dict[str, collections.Counter] = {}

    with tempfile.TemporaryDirectory() as td:
        for i, ((rnum, rname, _chain), aidx) in enumerate(pep_residues.items(), 1):
            tag = f"{rname}{rnum}"
            idx_expr = " ".join(str(a) for a in aidx)

            # --- contacts with the receptor -------------------------------
            size_xvg = os.path.join(td, f"n_{tag}.xvg")
            sel = (f'group "{rec_group}" and within {cutoff} of '
                   f'(atomnr {idx_expr})')
            try:
                run_gmx([gmx, "select", "-s", tpr, "-f", xtc, "-n", ndx,
                         "-select", sel, "-os", size_xvg] + time_args)
                rows = parse_xvg(size_xvg)
                counts = [r[1] for r in rows if len(r) > 1]
            except RuntimeError as exc:
                print(f"    WARNING: {tag}: {exc}".split("\n")[0])
                counts = []

            n_frames = len(counts)
            mean_c = sum(counts) / n_frames if n_frames else 0.0
            in_contact = sum(1 for c in counts if c > 0)
            inter_rows.append({
                "residue": tag,
                "resid": rnum,
                "resname": rname,
                "mean_contacts": mean_c,
                "max_contacts": int(max(counts)) if counts else 0,
                "frames_in_contact": in_contact,
                "occupancy": in_contact / n_frames if n_frames else 0.0,
                "n_frames": n_frames,
            })

            # --- which receptor residues does it touch? -------------------
            part_ndx = os.path.join(td, f"p_{tag}.ndx")
            try:
                run_gmx([gmx, "select", "-s", tpr, "-f", xtc, "-n", ndx,
                         "-select", sel, "-on", part_ndx,
                         "-selrpos", "atom"] + time_args)
                touched = collections.Counter()
                for gname, idxs in read_ndx_groups(part_ndx).items():
                    for a in idxs:
                        if a in atoms:
                            pr, pn, _, _pc = atoms[a]
                            touched[f"{pn}{pr}"] += 1
                partner_counts[tag] = touched
            except RuntimeError:
                partner_counts[tag] = collections.Counter()

            # --- intramolecular contacts (Fig 5B) -------------------------
            # exclude the residue itself and its sequence neighbours, which
            # are always "in contact" through the backbone
            excl = {rnum - 1, rnum, rnum + 1}
            others = [a for a in pep_atoms
                      if a in atoms and atoms[a][0] not in excl]
            if others:
                osize = os.path.join(td, f"i_{tag}.xvg")
                osel = (f'(atomnr {" ".join(str(a) for a in others)}) and '
                        f'within {cutoff} of (atomnr {idx_expr})')
                try:
                    run_gmx([gmx, "select", "-s", tpr, "-f", xtc, "-n", ndx,
                             "-select", osel, "-os", osize] + time_args)
                    orows = parse_xvg(osize)
                    ocounts = [r[1] for r in orows if len(r) > 1]
                except RuntimeError:
                    ocounts = []
                intra_rows.append({
                    "residue": tag,
                    "mean_contacts": (sum(ocounts) / len(ocounts)
                                      if ocounts else 0.0),
                    "max_contacts": int(max(ocounts)) if ocounts else 0,
                })
            else:
                intra_rows.append({"residue": tag, "mean_contacts": 0.0,
                                   "max_contacts": 0})

            print(f"    [{i:3d}/{len(pep_residues)}] {tag:<8s} "
                  f"inter {mean_c:7.1f}   intra "
                  f"{intra_rows[-1]['mean_contacts']:7.1f}")

    # ------------------------------------------------------------------ #
    #  write results
    # ------------------------------------------------------------------ #
    f_inter = os.path.join(outdir, "contacts_peptide_receptor.dat")
    with open(f_inter, "w") as fh:
        fh.write(f"# Contacts between peptide residues and the receptor\n")
        fh.write(f"# cutoff = {cutoff} nm ({cutoff*10:.1f} A), "
                 f"paper Fig 5A\n")
        fh.write("# residue  mean_contacts  max  frames_in_contact  "
                 "occupancy  n_frames\n")
        for r in inter_rows:
            fh.write(f"{r['residue']:<10s} {r['mean_contacts']:12.3f} "
                     f"{r['max_contacts']:6d} {r['frames_in_contact']:10d} "
                     f"{r['occupancy']:12.3f} {r['n_frames']:8d}\n")

    f_intra = os.path.join(outdir, "contacts_peptide_intra.dat")
    with open(f_intra, "w") as fh:
        fh.write("# Intramolecular contacts inside the peptide "
                 f"(cutoff {cutoff} nm), paper Fig 5B\n")
        fh.write("# i,i+-1 neighbours excluded\n")
        fh.write("# residue  mean_contacts  max\n")
        for r in intra_rows:
            fh.write(f"{r['residue']:<10s} {r['mean_contacts']:12.3f} "
                     f"{r['max_contacts']:6d}\n")

    f_part = os.path.join(outdir, "contact_partners.dat")
    with open(f_part, "w") as fh:
        fh.write("# Most frequent receptor partners for each peptide residue\n")
        fh.write("# (paper Table 1: contacts formed more than 10 times)\n")
        for tag, ctr in partner_counts.items():
            top = ctr.most_common(15)
            if not top:
                continue
            fh.write(f"\n{tag}:\n")
            for name, n in top:
                if n > 10:
                    fh.write(f"    {name:<10s} {n:8d}\n")

    # ranked summary, mirroring the paper's narrative
    ranked = sorted(inter_rows, key=lambda r: r["mean_contacts"], reverse=True)
    print()
    print("  " + "-" * 58)
    print("  Peptide residues ranked by contacts with the receptor")
    print("  " + "-" * 58)
    for r in ranked[:12]:
        print(f"    {r['residue']:<8s} {r['mean_contacts']:8.1f} contacts   "
              f"occupancy {r['occupancy']*100:5.1f}%")
    print("  " + "-" * 58)
    print(f"  wrote {f_inter}")
    print(f"  wrote {f_intra}")
    print(f"  wrote {f_part}")
    return 0


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Peptide-receptor and peptide-internal contact analysis")
    ap.add_argument("-s", "--tpr", required=True)
    ap.add_argument("-f", "--traj", required=True)
    ap.add_argument("-n", "--ndx", required=True)
    ap.add_argument("-o", "--outdir", default="analysis")
    ap.add_argument("--cutoff", type=float, default=0.7,
                    help="contact cut-off in nm (paper: 7 A = 0.7 nm)")
    ap.add_argument("--peptide-group", default="Peptide")
    ap.add_argument("--receptor-group", default="Receptor")
    ap.add_argument("-b", "--begin", type=float, default=None)
    ap.add_argument("-e", "--end", type=float, default=None)
    ap.add_argument("--gmx", default="gmx")
    args = ap.parse_args(argv)

    print("=" * 64)
    print("contacts.py  (paper Fig 5 / Table 1)")
    print(f"  cutoff {args.cutoff} nm = {args.cutoff*10:.1f} A")
    print("=" * 64)

    try:
        return analyse(args.gmx, args.tpr, args.traj, args.ndx, args.cutoff,
                       args.outdir, args.peptide_group, args.receptor_group,
                       args.begin, args.end)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
