#!/usr/bin/env bash
# =====================================================================
# analyze.sh - reproduce the analyses of the reference paper
#
#   Fig 1   backbone RMSD + RMSF of complex / receptor / peptide
#   Fig 2   RDF of the peptide around the receptor, whole run and in 4 quarters
#   Fig 3   SASA of the complex (LCPO in Amber -> gmx sasa here)
#   Fig 4   secondary-structure propensity of the peptide (helix/turn/bend)
#   Fig 5   contacts between peptide residues and the receptor (7 A cut-off)
#   Fig 6   water-mediated bridges  (scripts/water_bridges.py)
#   +       hydrogen bonds peptide<->receptor and inside the peptide
#   +       minimum distance to the PAS residues
#
# Everything lands in <workdir>/analysis/ as .xvg/.dat plus a summary.txt.
#
# Usage:
#   scripts/analyze.sh -d work/alllhrc
#   scripts/analyze.sh -d work/alllhrc --stage md --skip-bridges
# =====================================================================
set -euo pipefail

GMX="${GMX:-gmx}"
WORKDIR=""
STAGE="md"
OUTDIR=""
SKIP_BRIDGES=0
BRIDGE_DT=0
DT=""

usage() {
    sed -n '2,18p' "$0" | sed 's/^# \{0,1\}//'
    cat <<EOF

Options:
  -d, --dir DIR         system directory (required)
      --stage NAME      trajectory prefix to analyse (default: $STAGE)
  -o, --out DIR         output directory (default: <dir>/analysis)
      --skip-bridges    skip the water-bridge analysis (it is the slow one)
      --bridge-dt PS    only analyse frames every PS ps for bridges (0=all)
      --dt PS           only use frames every PS picoseconds
  -h, --help            this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dir)          WORKDIR="$2"; shift 2 ;;
        --stage)           STAGE="$2"; shift 2 ;;
        -o|--out)          OUTDIR="$2"; shift 2 ;;
        --skip-bridges)    SKIP_BRIDGES=1; shift ;;
        --bridge-dt)       BRIDGE_DT="$2"; shift 2 ;;
        --dt)              DT="$2"; shift 2 ;;
        -h|--help)         usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -n "$WORKDIR" ]] || { echo "ERROR: -d/--dir is required" >&2; usage; exit 1; }
[[ -d "$WORKDIR" ]] || { echo "ERROR: no such directory: $WORKDIR" >&2; exit 1; }

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKDIR="$(cd "$WORKDIR" && pwd)"
cd "$WORKDIR"

TPR="${STAGE}.tpr"
XTC="${STAGE}.xtc"
[[ -f "$TPR" ]] || { echo "ERROR: $TPR not found" >&2; exit 1; }
[[ -f "$XTC" ]] || { echo "ERROR: $XTC not found" >&2; exit 1; }
[[ -f index.ndx ]] || { echo "ERROR: index.ndx not found" >&2; exit 1; }

OUTDIR="${OUTDIR:-$WORKDIR/analysis}"
mkdir -p "$OUTDIR"

DT_ARGS=()
[[ -n "$DT" ]] && DT_ARGS+=(-dt "$DT")

echo "======================================================================"
echo " analyze.sh"
echo "   workdir : $WORKDIR"
echo "   traj    : $XTC   (tpr: $TPR)"
echo "   output  : $OUTDIR"
echo "======================================================================"

# ---------------------------------------------------------------------
# 0. make a PBC-corrected, receptor-fitted trajectory once and reuse it.
#    -pbc mol -center keeps the peptide next to the receptor instead of
#    jumping across the periodic boundary, which would ruin RMSD/contacts.
# ---------------------------------------------------------------------
FIT="${STAGE}_fit.xtc"
if [[ ! -f "$FIT" ]]; then
    echo
    echo ">>> [0] PBC correction + fit on the receptor"
    echo -e "Complex\nSystem" | "$GMX" trjconv -s "$TPR" -f "$XTC" -n index.ndx \
        -o "${STAGE}_whole.xtc" -pbc mol -center "${DT_ARGS[@]}"
    echo -e "Receptor_backbone\nSystem" | "$GMX" trjconv -s "$TPR" \
        -f "${STAGE}_whole.xtc" -n index.ndx -o "$FIT" -fit rot+trans
fi

# a matching reference structure (first frame, same treatment)
if [[ ! -f "${STAGE}_ref.pdb" ]]; then
    echo -e "Complex\nSystem" | "$GMX" trjconv -s "$TPR" -f "$XTC" -n index.ndx \
        -o "${STAGE}_ref.pdb" -pbc mol -center -dump 0
fi

# ---------------------------------------------------------------------
# 1. Fig 1 - backbone RMSD and RMSF
# ---------------------------------------------------------------------
echo
echo ">>> [1] RMSD / RMSF  (Fig 1)"
for grp in Complex_backbone Receptor_backbone Peptide_backbone; do
    echo -e "${grp}\n${grp}" | "$GMX" rms -s "$TPR" -f "$FIT" -n index.ndx \
        -o "$OUTDIR/rmsd_${grp}.xvg" -tu ns || true
done
for grp in Receptor_backbone Peptide_backbone; do
    echo "${grp}" | "$GMX" rmsf -s "$TPR" -f "$FIT" -n index.ndx \
        -o "$OUTDIR/rmsf_${grp}.xvg" -res || true
done

# ---------------------------------------------------------------------
# 2. Fig 2 - RDF of the peptide around the receptor, whole + 4 quarters
# ---------------------------------------------------------------------
echo
echo ">>> [2] RDF peptide around receptor  (Fig 2)"
echo -e "Receptor\nPeptide" | "$GMX" rdf -s "$TPR" -f "$FIT" -n index.ndx \
    -o "$OUTDIR/rdf_peptide_around_receptor.xvg" -bin 0.02 -rmax 3.0 || true

# split into four equal quarters (paper: "we divided the production dynamics
# trajectory into four equal parts and estimated the RDF for each of them")
TEND=$("$GMX" check -f "$FIT" 2>&1 | tr '\r' '\n' \
       | awk '/Last frame/{print $NF}' | tail -1)
if [[ -n "${TEND:-}" ]]; then
    for i in 1 2 3 4; do
        B=$(awk -v t="$TEND" -v i="$i" 'BEGIN{printf "%.4f", (i-1)*t/4.0}')
        E=$(awk -v t="$TEND" -v i="$i" 'BEGIN{printf "%.4f", i*t/4.0}')
        echo "    quarter $i: $B - $E ps"
        echo -e "Receptor\nPeptide" | "$GMX" rdf -s "$TPR" -f "$FIT" -n index.ndx \
            -o "$OUTDIR/rdf_quarter${i}.xvg" -bin 0.02 -rmax 3.0 \
            -b "$B" -e "$E" || true
    done
else
    echo "    WARNING: could not determine the trajectory end time; "
    echo "             skipping the four-quarter RDF split"
fi

# ---------------------------------------------------------------------
# 3. Fig 3 - SASA of the complex
# ---------------------------------------------------------------------
echo
echo ">>> [3] SASA  (Fig 3)"
"$GMX" sasa -s "$TPR" -f "$FIT" -n index.ndx \
    -surface 'group "Complex"' -output 'group "Peptide"' \
    -o "$OUTDIR/sasa_complex.xvg" -tu ns || true

# ---------------------------------------------------------------------
# 4. Fig 4 - secondary structure of the peptide
# ---------------------------------------------------------------------
echo
echo ">>> [4] secondary structure of the peptide  (Fig 4)"
if "$GMX" dssp -h >/dev/null 2>&1; then
    "$GMX" dssp -s "$TPR" -f "$FIT" -n index.ndx -sel 'group "Peptide"' \
        -o "$OUTDIR/dssp_peptide.dat" -num "$OUTDIR/dssp_count.xvg" || true
elif "$GMX" do_dssp -h >/dev/null 2>&1; then
    echo "Peptide" | "$GMX" do_dssp -s "$TPR" -f "$FIT" -n index.ndx \
        -o "$OUTDIR/dssp_peptide.xpm" -sc "$OUTDIR/dssp_count.xvg" || true
else
    echo "    (no dssp module in this GROMACS build - skipped)"
fi

# ---------------------------------------------------------------------
# 5. Fig 5 - contacts within 7 A, per peptide residue
# ---------------------------------------------------------------------
echo
echo ">>> [5] contacts peptide<->receptor, 7 A  (Fig 5)"
python3 "$REPO/scripts/contacts.py" -s "$TPR" -f "$FIT" -n index.ndx \
    --cutoff 0.7 --gmx "$GMX" -o "$OUTDIR" || true

# ---------------------------------------------------------------------
# 6. hydrogen bonds
# ---------------------------------------------------------------------
echo
echo ">>> [6] hydrogen bonds"
echo -e "Peptide\nReceptor" | "$GMX" hbond -s "$TPR" -f "$FIT" -n index.ndx \
    -num "$OUTDIR/hbond_peptide_receptor.xvg" || true
echo -e "Peptide\nPeptide" | "$GMX" hbond -s "$TPR" -f "$FIT" -n index.ndx \
    -num "$OUTDIR/hbond_peptide_intra.xvg" || true

# ---------------------------------------------------------------------
# 7. distance to the PAS
# ---------------------------------------------------------------------
echo
echo ">>> [7] minimum distance peptide <-> PAS"
if grep -q "^\[ PAS \]" index.ndx 2>/dev/null; then
    echo -e "Peptide\nPAS" | "$GMX" mindist -s "$TPR" -f "$FIT" -n index.ndx \
        -od "$OUTDIR/mindist_peptide_PAS.xvg" -tu ns || true
else
    echo "    (no PAS group in index.ndx - add one with scripts/make_pas_index.sh)"
fi
echo -e "Peptide\nReceptor" | "$GMX" mindist -s "$TPR" -f "$FIT" -n index.ndx \
    -od "$OUTDIR/mindist_peptide_receptor.xvg" -tu ns || true

# ---------------------------------------------------------------------
# 8. Fig 6 - water-mediated bridges
# ---------------------------------------------------------------------
if [[ $SKIP_BRIDGES -eq 0 ]]; then
    echo
    echo ">>> [8] water-mediated bridges  (Fig 6)"
    python3 "$REPO/scripts/water_bridges.py" -s "$TPR" -f "$FIT" -n index.ndx \
        --dt "$BRIDGE_DT" --gmx "$GMX" -o "$OUTDIR" || true
else
    echo
    echo ">>> [8] water bridges SKIPPED (--skip-bridges)"
fi

# ---------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------
{
    echo "Analysis summary - $(date)"
    echo "workdir: $WORKDIR"
    echo "traj   : $XTC"
    echo
    for f in "$OUTDIR"/*.xvg; do
        [[ -f "$f" ]] || continue
        n=$(grep -vc '^[#@]' "$f" || echo 0)
        printf "  %-45s %6s data rows\n" "$(basename "$f")" "$n"
    done
} > "$OUTDIR/summary.txt"

echo
echo "======================================================================"
echo " ANALYSIS DONE -> $OUTDIR"
cat "$OUTDIR/summary.txt" | tail -n +5
echo "======================================================================"
