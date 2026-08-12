#!/usr/bin/env bash
# =====================================================================
# build_system.sh - system setup for the AChE(4EY6)-peptide PAS complex
#
# Reproduces section 2.2 of the reference paper:
#   * caps on every chain end and every chain break
#   * protonation at pH 7
#   * TIP3P water in a TRUNCATED OCTAHEDRON box
#   * NaCl at physiological concentration + neutralisation
#
# Produces, in <workdir>:
#   prepared.pdb  cleaned + capped complex        (prepare_structure.py)
#   topol.top     ff14SB topology
#   solv_ions.gro coordinates ready for minimisation
#   index.ndx     groups: Protein_Peptide / Receptor / Peptide / Water_and_ions
#
# Usage:
#   scripts/build_system.sh -f input/alllhrc_complex.pdb -d work/alllhrc
#   scripts/build_system.sh -f in.pdb -d work/x --peptide-chain B --conc 0.15
# =====================================================================
set -euo pipefail

# ---------------------------------------------------------------- defaults
GMX="${GMX:-gmx}"
INPUT=""
WORKDIR=""
PEPTIDE_CHAIN=""
BOXDIST="1.2"          # nm, min solute-box distance
CONC="0.15"            # mol/L NaCl, physiological
FF="amber14sb"
WATER="tip3p"
BOXTYPE="octahedron"   # paper: truncated octahedron
FREE_TERMINI=""
EXTRA_PREP_ARGS=""

usage() {
    sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'
    cat <<EOF

Options:
  -f, --file PDB          input docked complex (required)
  -d, --dir DIR           working directory (required)
      --peptide-chain ID  chain ID of the peptide (default: last chain seen)
      --box-dist NM       solute-to-box distance      (default: $BOXDIST)
      --conc M            NaCl concentration in mol/L (default: $CONC)
      --box-type TYPE     octahedron | cubic | dodecahedron (default: $BOXTYPE)
      --free-termini IDS  comma-separated chains whose real termini stay charged
      --ff NAME           force field directory name  (default: $FF)
      --water NAME        water model                 (default: $WATER)
  -h, --help              this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--file)          INPUT="$2"; shift 2 ;;
        -d|--dir)           WORKDIR="$2"; shift 2 ;;
        --peptide-chain)    PEPTIDE_CHAIN="$2"; shift 2 ;;
        --box-dist)         BOXDIST="$2"; shift 2 ;;
        --conc)             CONC="$2"; shift 2 ;;
        --box-type)         BOXTYPE="$2"; shift 2 ;;
        --free-termini)     FREE_TERMINI="$2"; shift 2 ;;
        --ff)               FF="$2"; shift 2 ;;
        --water)            WATER="$2"; shift 2 ;;
        -h|--help)          usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -n "$INPUT"   ]] || { echo "ERROR: -f/--file is required" >&2; usage; exit 1; }
[[ -n "$WORKDIR" ]] || { echo "ERROR: -d/--dir is required"  >&2; usage; exit 1; }
[[ -f "$INPUT"   ]] || { echo "ERROR: no such file: $INPUT"  >&2; exit 1; }

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT="$(cd "$(dirname "$INPUT")" && pwd)/$(basename "$INPUT")"

command -v "$GMX" >/dev/null 2>&1 || {
    echo "ERROR: '$GMX' not found. source your GMXRC or set GMX=..." >&2; exit 1; }

mkdir -p "$WORKDIR"
WORKDIR="$(cd "$WORKDIR" && pwd)"
cd "$WORKDIR"

echo "======================================================================"
echo " build_system.sh"
echo "   input     : $INPUT"
echo "   workdir   : $WORKDIR"
echo "   force fld : $FF   water: $WATER"
echo "   box       : $BOXTYPE, $BOXDIST nm   NaCl: $CONC M"
echo "======================================================================"

# pdb2gmx searches the CWD for <name>.ff before $GMXDATA/top
if [[ -d "$REPO/forcefield/${FF}.ff" && ! -e "${FF}.ff" ]]; then
    cp -r "$REPO/forcefield/${FF}.ff" .
    echo ">>> vendored force field ${FF}.ff copied into workdir"
fi

# ---------------------------------------------------------------- 1. clean
echo
echo ">>> [1/6] cleaning + capping (ACE/NME on all chain ends and breaks)"
PREP_ARGS=(-o prepared.pdb)
[[ -n "$PEPTIDE_CHAIN" ]] && PREP_ARGS+=(--peptide-chain "$PEPTIDE_CHAIN")
[[ -n "$FREE_TERMINI"  ]] && PREP_ARGS+=(--free-termini  "$FREE_TERMINI")
python3 "$REPO/scripts/prepare_structure.py" "$INPUT" "${PREP_ARGS[@]}"

N_SEG=$(python3 -c "import json;print(json.load(open('prepared.json'))['n_segments'])")
N_TER=$(python3 -c "import json;print(json.load(open('prepared.json'))['pdb2gmx_ter_answers'])")

# ---------------------------------------------------------------- 2. topology
echo
echo ">>> [2/6] pdb2gmx  ($N_SEG segments -> $N_TER terminus prompts, all 'None')"
# Every terminus is a named cap (ACE/NME), so the answer is always the LAST
# option in the list, which pdb2gmx labels "None".  Feed it that many times.
# (built with a loop rather than `yes | head` so no SIGPIPE trips pipefail)
TER_ANSWERS=""
for _ in $(seq "$N_TER"); do TER_ANSWERS+="None"$'\n'; done
printf '%s' "$TER_ANSWERS" | \
    "$GMX" pdb2gmx -f prepared.pdb -o processed.gro -p topol.top -i posre.itp \
        -ff "$FF" -water "$WATER" -ter -ignh -merge no

# paper: 3 kcal/mol/A^2 harmonic restraints on protein+peptide heavy atoms
echo
echo ">>> [2b/6] setting restraints to 3 kcal/mol/A^2 (-DPOSRES_HEAVY)"
python3 "$REPO/scripts/make_posre.py" --dir "$WORKDIR"

# ---------------------------------------------------------------- 3. box
echo
echo ">>> [3/6] editconf: $BOXTYPE box, $BOXDIST nm from solute"
case "$BOXTYPE" in
    octahedron)   BT_ARGS=(-bt octahedron) ;;
    cubic)        BT_ARGS=(-bt cubic) ;;
    dodecahedron) BT_ARGS=(-bt dodecahedron) ;;
    *) echo "ERROR: unknown box type $BOXTYPE" >&2; exit 1 ;;
esac
"$GMX" editconf -f processed.gro -o boxed.gro -c -d "$BOXDIST" "${BT_ARGS[@]}"

# ---------------------------------------------------------------- 4. solvate
echo
echo ">>> [4/6] solvate with TIP3P"
"$GMX" solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

# ---------------------------------------------------------------- 5. ions
echo
echo ">>> [5/6] adding NaCl to $CONC M + neutralising"
cat > ions.mdp <<'EOF'
; minimal mdp just to produce a .tpr for genion
integrator      = steep
nsteps          = 0
cutoff-scheme   = Verlet
coulombtype     = PME
rcoulomb        = 1.2
rvdw            = 1.2
pbc             = xyz
EOF
"$GMX" grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 5

# "SOL" is the group genion replaces
echo SOL | "$GMX" genion -s ions.tpr -o solv_ions.gro -p topol.top \
    -pname NA -nname CL -neutral -conc "$CONC"

# ---------------------------------------------------------------- 6. index
echo
echo ">>> [6/6] building index groups"
python3 "$REPO/scripts/make_index.py" \
    --struct solv_ions.gro --top topol.top --prepared-json prepared.json \
    --gmx "$GMX" -o index.ndx

# ---------------------------------------------------------------- summary
NATOMS=$(head -2 solv_ions.gro | tail -1 | tr -d ' ')
echo
echo "======================================================================"
echo " SYSTEM READY in $WORKDIR"
echo "   total atoms : $NATOMS      (paper: >76,000)"
echo "   topology    : topol.top"
echo "   coordinates : solv_ions.gro"
echo "   index       : index.ndx"
echo
echo " Next:  scripts/run_md.sh -d $WORKDIR            # production"
echo "        scripts/run_md.sh -d $WORKDIR --test     # 100-step smoke test"
echo "======================================================================"
