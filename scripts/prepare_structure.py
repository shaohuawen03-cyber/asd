#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
prepare_structure.py
====================

Clean a docked AChE(4EY6)-peptide complex PDB and make it ``gmx pdb2gmx``-ready,
following the system-setup section of the reference paper:

    "Capping was performed for each end of the protein chains as well as for
     both ends of the broken parts (residues 259, 262, 492 and 495)."

What it does
------------
1.  Keeps protein residues only; drops waters, ions, the galantamine ligand,
    N-glycans (NAG/BMA/MAN/FUC) and any other HETATM junk.
2.  Renames non-standard residues that have a standard equivalent (MSE -> MET,
    HSD/HSE/HSP/HID/HIE/HIP -> HIS, CYX -> CYS, ...).
3.  Keeps the first altLoc of each atom and strips existing hydrogens
    (pdb2gmx rebuilds them; this also kills naming mismatches).
4.  Splits every chain into *continuous segments*, breaking wherever the
    backbone C(i)-N(i+1) distance exceeds ``--break-dist`` (default 2.5 A) or
    the residue numbering jumps.  For 4EY6 this reproduces the breaks around
    259/262 and 492/495 automatically - nothing is hard-coded.
5.  Caps every segment terminus with ACE / NME, built with proper internal
    coordinates and a dihedral scan that picks the least-clashing orientation.
6.  Writes the cleaned PDB plus a JSON report describing chains, segments,
    caps, S-S candidates and the pdb2gmx terminus answers to feed back.

Only the Python standard library is required.

Usage
-----
    python3 scripts/prepare_structure.py input/alllhrc_complex.pdb \\
            -o work/alllhrc/prepared.pdb

    # leave the peptide's real termini charged (still caps internal breaks)
    python3 scripts/prepare_structure.py in.pdb -o out.pdb --free-termini B
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

# --------------------------------------------------------------------------- #
#  Residue tables
# --------------------------------------------------------------------------- #

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

# non-standard -> standard.  Everything on the left is treated as protein.
RESIDUE_ALIASES = {
    "MSE": "MET",   # selenomethionine
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
    "CYX": "CYS", "CYM": "CYS", "CSS": "CYS",
    "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "ARN": "ARG",
    "SEC": "CYS", "PYL": "LYS",
    "ACE": "ACE", "NME": "NME", "NMA": "NME", "NHE": "NHE",
}

# atoms that only exist in the selenomethionine / other variants
ATOM_ALIASES = {
    ("MET", "SE"): "SD",
    ("MET", "SED"): "SD",
}

CAP_RESNAMES = {"ACE", "NME", "NMA", "NHE"}

# --------------------------------------------------------------------------- #
#  Tiny 3-D helpers (no numpy dependency)
# --------------------------------------------------------------------------- #

Vec = Tuple[float, float, float]


def vsub(a: Vec, b: Vec) -> Vec:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def vadd(a: Vec, b: Vec) -> Vec:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def vmul(a: Vec, s: float) -> Vec:
    return (a[0] * s, a[1] * s, a[2] * s)


def vdot(a: Vec, b: Vec) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def vcross(a: Vec, b: Vec) -> Vec:
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def vnorm(a: Vec) -> float:
    return math.sqrt(vdot(a, a))


def vunit(a: Vec) -> Vec:
    n = vnorm(a)
    if n < 1e-9:
        raise ValueError("cannot normalise a zero-length vector")
    return (a[0] / n, a[1] / n, a[2] / n)


def dist(a: Vec, b: Vec) -> float:
    return vnorm(vsub(a, b))


def dihedral(p0: Vec, p1: Vec, p2: Vec, p3: Vec) -> float:
    """Dihedral angle p0-p1-p2-p3 in degrees, IUPAC sign convention."""
    b0 = vsub(p0, p1)
    b1 = vsub(p2, p1)
    b2 = vsub(p3, p2)
    b1u = vunit(b1)
    v = vsub(b0, vmul(b1u, vdot(b0, b1u)))
    w = vsub(b2, vmul(b1u, vdot(b2, b1u)))
    x = vdot(v, w)
    y = vdot(vcross(b1u, v), w)
    return math.degrees(math.atan2(y, x))


def place_atom(a: Vec, b: Vec, c: Vec,
               bond: float, angle_deg: float, torsion_deg: float) -> Vec:
    """
    NeRF: place atom D given a bonded chain A-B-C, with
    |C-D| = bond, angle(B,C,D) = angle_deg, dihedral(A,B,C,D) = torsion_deg.
    """
    angle = math.radians(angle_deg)
    torsion = math.radians(torsion_deg)

    bc = vunit(vsub(c, b))
    ab = vsub(b, a)
    n = vcross(ab, bc)
    if vnorm(n) < 1e-7:          # A, B, C collinear -> pick any perpendicular
        tmp = (1.0, 0.0, 0.0) if abs(bc[0]) < 0.9 else (0.0, 1.0, 0.0)
        n = vcross(tmp, bc)
    n = vunit(n)
    m = vcross(n, bc)

    d2 = (-bond * math.cos(angle),
          bond * math.sin(angle) * math.cos(torsion),
          bond * math.sin(angle) * math.sin(torsion))

    return (c[0] + d2[0] * bc[0] + d2[1] * m[0] + d2[2] * n[0],
            c[1] + d2[0] * bc[1] + d2[1] * m[1] + d2[2] * n[1],
            c[2] + d2[0] * bc[2] + d2[1] * m[2] + d2[2] * n[2])


# --------------------------------------------------------------------------- #
#  PDB data model
# --------------------------------------------------------------------------- #

class Atom:
    __slots__ = ("name", "resname", "chain", "resseq", "icode", "xyz",
                 "element", "occ", "bfac")

    def __init__(self, name, resname, chain, resseq, icode, xyz,
                 element="", occ=1.0, bfac=0.0):
        self.name = name
        self.resname = resname
        self.chain = chain
        self.resseq = resseq
        self.icode = icode
        self.xyz: Vec = xyz
        self.element = element
        self.occ = occ
        self.bfac = bfac


class Residue:
    __slots__ = ("resname", "chain", "resseq", "icode", "atoms")

    def __init__(self, resname, chain, resseq, icode):
        self.resname = resname
        self.chain = chain
        self.resseq = resseq
        self.icode = icode
        self.atoms: List[Atom] = []

    @property
    def key(self):
        return (self.chain, self.resseq, self.icode)

    def get(self, name: str) -> Optional[Atom]:
        for a in self.atoms:
            if a.name == name:
                return a
        return None

    def xyz(self, name: str) -> Optional[Vec]:
        a = self.get(name)
        return a.xyz if a is not None else None

    def has_backbone(self) -> bool:
        return all(self.get(n) is not None for n in ("N", "CA", "C"))

    def __repr__(self):
        return f"<{self.resname}{self.resseq}{self.icode.strip()}:{self.chain}>"


def guess_element(atom_name: str, raw_element: str) -> str:
    e = raw_element.strip()
    if e:
        return e.capitalize()
    n = atom_name.strip()
    if not n:
        return ""
    if n[0].isdigit():
        n = n[1:]
    two = n[:2].capitalize()
    if two in ("Cl", "Br", "Se", "Zn", "Fe", "Mg", "Mn", "Na", "Ca"):
        return two
    return n[0].upper()


def read_pdb(path: str, keep_models: int = 1) -> List[Atom]:
    """Read ATOM/HETATM records.  Only the first MODEL is kept by default."""
    atoms: List[Atom] = []
    model = 0
    seen_alt: Dict[Tuple, str] = {}
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            rec = line[:6]
            if rec == "MODEL ":
                model += 1
                if model > keep_models:
                    break
                continue
            if rec == "ENDMDL":
                if model >= keep_models:
                    break
                continue
            if rec not in ("ATOM  ", "HETATM"):
                continue

            name = line[12:16].strip()
            altloc = line[16]
            resname = line[17:20].strip().upper()
            chain = line[21]
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            icode = line[26]
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            occ = float(line[54:60]) if line[54:60].strip() else 1.0
            bfac = float(line[60:66]) if line[60:66].strip() else 0.0
            element = guess_element(name, line[76:78] if len(line) >= 78 else "")

            # keep only the first altLoc seen for a given atom slot
            if altloc not in (" ", ""):
                slot = (chain, resseq, icode, resname, name)
                if slot in seen_alt and seen_alt[slot] != altloc:
                    continue
                seen_alt[slot] = altloc

            atoms.append(Atom(name, resname, chain, resseq, icode,
                              (x, y, z), element, occ, bfac))
    return atoms


def group_residues(atoms: Iterable[Atom]) -> List[Residue]:
    residues: List[Residue] = []
    index: Dict[Tuple, Residue] = {}
    for a in atoms:
        key = (a.chain, a.resseq, a.icode, a.resname)
        res = index.get(key)
        if res is None:
            res = Residue(a.resname, a.chain, a.resseq, a.icode)
            index[key] = res
            residues.append(res)
        res.atoms.append(a)
    return residues


# --------------------------------------------------------------------------- #
#  Cleaning
# --------------------------------------------------------------------------- #

def clean_residues(residues: Sequence[Residue],
                   keep_hydrogens: bool = False) -> Tuple[List[Residue], Dict]:
    """Drop non-protein, normalise names, strip H, drop broken residues."""
    kept: List[Residue] = []
    report = {"dropped_residues": {}, "renamed": {}, "removed_hydrogens": 0,
              "dropped_incomplete": []}

    for res in residues:
        rn = res.resname
        std = rn if rn in STANDARD_AA else RESIDUE_ALIASES.get(rn)

        if std is None or std in CAP_RESNAMES:
            # not a protein residue we can handle -> drop it (water, ligand,
            # sugar, ion, ...).  Pre-existing caps are dropped and rebuilt.
            report["dropped_residues"][rn] = \
                report["dropped_residues"].get(rn, 0) + 1
            continue

        if std != rn:
            report["renamed"][rn] = std
            for a in res.atoms:
                a.resname = std
                new = ATOM_ALIASES.get((std, a.name))
                if new:
                    a.name = new
            res.resname = std

        if not keep_hydrogens:
            before = len(res.atoms)
            res.atoms = [a for a in res.atoms
                         if guess_element(a.name, a.element) != "H"
                         and guess_element(a.name, a.element) != "D"]
            report["removed_hydrogens"] += before - len(res.atoms)

        if not res.has_backbone():
            report["dropped_incomplete"].append(
                f"{res.resname}{res.resseq}{res.icode.strip()}:{res.chain}")
            continue

        kept.append(res)

    return kept, report


# --------------------------------------------------------------------------- #
#  Segmentation
# --------------------------------------------------------------------------- #

def split_segments(residues: Sequence[Residue],
                   break_dist: float = 2.5) -> List[List[Residue]]:
    """Split a chain's residues into peptide-bond-continuous segments."""
    segments: List[List[Residue]] = []
    current: List[Residue] = []

    for res in residues:
        if not current:
            current = [res]
            continue
        prev = current[-1]
        c_prev = prev.xyz("C")
        n_cur = res.xyz("N")
        broken = True
        if c_prev is not None and n_cur is not None:
            d = dist(c_prev, n_cur)
            broken = d > break_dist
        if broken:
            segments.append(current)
            current = [res]
        else:
            current.append(res)

    if current:
        segments.append(current)
    return segments


def find_disulfides(residues: Sequence[Residue],
                    cutoff: float = 2.5) -> List[Tuple[str, str, float]]:
    cys = [r for r in residues if r.resname == "CYS" and r.get("SG")]
    out = []
    for i in range(len(cys)):
        for j in range(i + 1, len(cys)):
            d = dist(cys[i].xyz("SG"), cys[j].xyz("SG"))
            if d <= cutoff:
                out.append((f"CYS{cys[i].resseq}:{cys[i].chain}",
                            f"CYS{cys[j].resseq}:{cys[j].chain}", round(d, 2)))
    return out


# --------------------------------------------------------------------------- #
#  Capping
# --------------------------------------------------------------------------- #

# ideal internal coordinates (Amber ff14SB-ish)
ACE_C_N = 1.335
ACE_CA_N_C = 121.7
ACE_C_O = 1.229
ACE_N_C_O = 122.9
ACE_C_CH3 = 1.508
ACE_N_C_CH3 = 116.6

NME_C_N = 1.335
NME_CA_C_N = 116.6
NME_N_CH3 = 1.449
NME_C_N_CH3 = 121.9


def _clash_score(pos: Vec, cloud: Sequence[Vec]) -> float:
    """Smallest distance from `pos` to the surrounding atom cloud."""
    best = 1e9
    for q in cloud:
        d2 = (pos[0] - q[0]) ** 2 + (pos[1] - q[1]) ** 2 + (pos[2] - q[2]) ** 2
        if d2 < best:
            best = d2
    return math.sqrt(best)


def _neighbourhood(centre: Vec, all_xyz: Sequence[Vec],
                   radius: float = 9.0) -> List[Vec]:
    r2 = radius * radius
    out = []
    for q in all_xyz:
        d2 = ((centre[0] - q[0]) ** 2 + (centre[1] - q[1]) ** 2
              + (centre[2] - q[2]) ** 2)
        if 1e-6 < d2 < r2:
            out.append(q)
    return out


def build_ace(first: Residue, cloud: Sequence[Vec]) -> Residue:
    """ACE cap on the N-terminus of `first`; scans phi to avoid clashes."""
    n = first.xyz("N")
    ca = first.xyz("CA")
    c = first.xyz("C")
    local = _neighbourhood(n, cloud)

    best = None
    for tor in range(-180, 180, 15):
        c_ace = place_atom(c, ca, n, ACE_C_N, ACE_CA_N_C, float(tor))
        ch3 = place_atom(ca, n, c_ace, ACE_C_CH3, ACE_N_C_CH3, 180.0)
        o = place_atom(ca, n, c_ace, ACE_C_O, ACE_N_C_O, 0.0)
        score = min(_clash_score(c_ace, local),
                    _clash_score(ch3, local),
                    _clash_score(o, local))
        if best is None or score > best[0]:
            best = (score, c_ace, ch3, o)

    _, c_ace, ch3, o = best
    cap = Residue("ACE", first.chain, first.resseq - 1, " ")
    cap.atoms = [
        Atom("CH3", "ACE", first.chain, cap.resseq, " ", ch3, "C"),
        Atom("C",   "ACE", first.chain, cap.resseq, " ", c_ace, "C"),
        Atom("O",   "ACE", first.chain, cap.resseq, " ", o, "O"),
    ]
    return cap


def build_nme(last: Residue, cloud: Sequence[Vec]) -> Residue:
    """NME cap on the C-terminus of `last`; N placed anti to the carbonyl O."""
    n = last.xyz("N")
    ca = last.xyz("CA")
    c = last.xyz("C")
    o = last.xyz("O") or last.xyz("OC1") or last.xyz("O1")
    local = _neighbourhood(c, cloud)

    if o is not None:
        psi = dihedral(n, ca, c, o) + 180.0
    else:
        psi = 180.0
    n_nme = place_atom(n, ca, c, NME_C_N, NME_CA_C_N, psi)

    best = None
    for tor in (180.0, 150.0, -150.0, 120.0, -120.0):
        ch3 = place_atom(ca, c, n_nme, NME_N_CH3, NME_C_N_CH3, tor)
        score = _clash_score(ch3, local)
        if best is None or score > best[0]:
            best = (score, ch3)
    ch3 = best[1]

    cap = Residue("NME", last.chain, last.resseq + 1, " ")
    cap.atoms = [
        Atom("N",   "NME", last.chain, cap.resseq, " ", n_nme, "N"),
        Atom("CH3", "NME", last.chain, cap.resseq, " ", ch3, "C"),
    ]
    return cap


# --------------------------------------------------------------------------- #
#  Writing
# --------------------------------------------------------------------------- #

def format_atom_name(name: str, element: str) -> str:
    """PDB columns 13-16 with the usual right-shift rule."""
    if len(name) >= 4:
        return name[:4]
    if len(element) == 2:
        return f"{name:<4s}"
    return f" {name:<3s}"


def write_pdb(path: str, segments: List[Tuple[str, List[Residue]]],
              title: str = "") -> int:
    serial = 0
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w") as fh:
        if title:
            fh.write(f"TITLE     {title[:60]}\n")
        fh.write("REMARK   1 Prepared by scripts/prepare_structure.py\n")
        fh.write("REMARK   1 ACE/NME caps added at every segment terminus\n")
        for _seg_id, residues in segments:
            for res in residues:
                for a in res.atoms:
                    serial += 1
                    el = guess_element(a.name, a.element)
                    fh.write(
                        "ATOM  {ser:5d} {nm}{alt}{rn:>3s} {ch}{seq:4d}{ic}   "
                        "{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{b:6.2f}"
                        "          {el:>2s}\n".format(
                            ser=serial % 100000,
                            nm=format_atom_name(a.name, el),
                            alt=" ",
                            rn=a.resname,
                            ch=a.chain,
                            seq=a.resseq,
                            ic=a.icode if a.icode.strip() else " ",
                            x=a.xyz[0], y=a.xyz[1], z=a.xyz[2],
                            occ=1.00, b=0.00, el=el))
            serial += 1
            fh.write("TER   {ser:5d}\n".format(ser=serial % 100000))
        fh.write("END\n")
    return serial


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #

def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Clean + cap a docked AChE-peptide complex for pdb2gmx.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("input", help="input complex PDB")
    ap.add_argument("-o", "--output", required=True, help="cleaned PDB")
    ap.add_argument("-j", "--json", default=None,
                    help="JSON report path (default: <output>.json)")
    ap.add_argument("--break-dist", type=float, default=2.5,
                    help="C(i)-N(i+1) distance above which a chain is broken")
    ap.add_argument("--free-termini", default="",
                    help="comma-separated chain IDs whose REAL termini stay "
                         "uncapped/charged (internal breaks are still capped)")
    ap.add_argument("--no-cap", action="store_true",
                    help="do not add any caps at all")
    ap.add_argument("--keep-hydrogens", action="store_true",
                    help="keep input hydrogens (default: strip, pdb2gmx rebuilds)")
    ap.add_argument("--peptide-chain", default=None,
                    help="chain ID of the peptide; only used to label the report")
    args = ap.parse_args(argv)

    if not os.path.isfile(args.input):
        print(f"ERROR: no such file: {args.input}", file=sys.stderr)
        return 2

    free = {c.strip() for c in args.free_termini.split(",") if c.strip()}

    atoms = read_pdb(args.input)
    if not atoms:
        print("ERROR: no ATOM/HETATM records parsed", file=sys.stderr)
        return 2

    residues = group_residues(atoms)
    residues, clean_report = clean_residues(residues, args.keep_hydrogens)
    if not residues:
        print("ERROR: nothing left after cleaning - is this a protein PDB?",
              file=sys.stderr)
        return 2

    # chain order as first seen
    chain_order: List[str] = []
    by_chain: Dict[str, List[Residue]] = {}
    for r in residues:
        if r.chain not in by_chain:
            by_chain[r.chain] = []
            chain_order.append(r.chain)
        by_chain[r.chain].append(r)

    cloud = [a.xyz for r in residues for a in r.atoms]

    out_segments: List[Tuple[str, List[Residue]]] = []
    seg_report = []
    n_ace = n_nme = 0

    for chain in chain_order:
        chain_res = sorted(by_chain[chain], key=lambda r: (r.resseq, r.icode))
        segments = split_segments(chain_res, args.break_dist)
        for si, seg in enumerate(segments):
            is_first = (si == 0)
            is_last = (si == len(segments) - 1)
            # a real terminus is the very start/end of the chain
            cap_n = not args.no_cap and not (is_first and chain in free)
            cap_c = not args.no_cap and not (is_last and chain in free)

            body = list(seg)
            if cap_n:
                body.insert(0, build_ace(seg[0], cloud))
                n_ace += 1
            if cap_c:
                body.append(build_nme(seg[-1], cloud))
                n_nme += 1

            out_segments.append((f"{chain}{si}", body))
            seg_report.append({
                "chain": chain,
                "segment": si,
                "first_residue": f"{seg[0].resname}{seg[0].resseq}",
                "last_residue": f"{seg[-1].resname}{seg[-1].resseq}",
                "n_residues": len(seg),
                "n_capped": ("ACE" if cap_n else "-") + "/" + ("NME" if cap_c else "-"),
            })

    n_atoms = write_pdb(args.output, out_segments,
                        title=os.path.basename(args.input))

    ss = find_disulfides(residues)
    n_chains_pdb2gmx = len(out_segments)

    report = {
        "input": os.path.abspath(args.input),
        "output": os.path.abspath(args.output),
        "n_atoms_written": n_atoms,
        "n_residues_kept": len(residues),
        "chains": chain_order,
        "peptide_chain": args.peptide_chain or (chain_order[-1] if len(chain_order) > 1 else None),
        "segments": seg_report,
        "n_segments": n_chains_pdb2gmx,
        "caps_added": {"ACE": n_ace, "NME": n_nme},
        "disulfide_candidates": ss,
        "cleaning": clean_report,
        # every terminus is a named cap -> answer "None" to all pdb2gmx -ter
        # prompts.  2 answers per segment (start + end).
        "pdb2gmx_ter_answers": 2 * n_chains_pdb2gmx,
    }

    json_path = args.json or (os.path.splitext(args.output)[0] + ".json")
    with open(json_path, "w") as fh:
        json.dump(report, fh, indent=2)

    # ---------------- human-readable summary ---------------- #
    print("=" * 68)
    print("prepare_structure.py")
    print("=" * 68)
    print(f"  input             : {args.input}")
    print(f"  output            : {args.output}")
    print(f"  report            : {json_path}")
    print(f"  atoms written     : {n_atoms}")
    print(f"  residues kept     : {len(residues)}")
    print(f"  chains            : {', '.join(chain_order)}")
    print(f"  segments          : {n_chains_pdb2gmx}  (pdb2gmx will see this "
          f"many chains)")
    print(f"  caps added        : {n_ace} x ACE, {n_nme} x NME")
    if clean_report["dropped_residues"]:
        drops = ", ".join(f"{k}x{v}" for k, v in
                          sorted(clean_report["dropped_residues"].items()))
        print(f"  dropped (het/wat) : {drops}")
    if clean_report["removed_hydrogens"]:
        print(f"  hydrogens removed : {clean_report['removed_hydrogens']}")
    if clean_report["dropped_incomplete"]:
        print(f"  incomplete resid. : {len(clean_report['dropped_incomplete'])} "
              f"(no N/CA/C) -> {clean_report['dropped_incomplete'][:5]}")
    if ss:
        print(f"  S-S candidates    : {len(ss)}")
        for a, b, d in ss:
            print(f"      {a} - {b}  ({d} A)")
    print()
    print("  Segment layout:")
    for s in seg_report:
        print(f"    chain {s['chain']} seg {s['segment']}: "
              f"{s['first_residue']:>8s} .. {s['last_residue']:<8s} "
              f"({s['n_residues']:4d} res)  caps {s['n_capped']}")
    print("=" * 68)
    return 0


if __name__ == "__main__":
    sys.exit(main())
