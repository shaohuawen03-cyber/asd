#!/usr/bin/env bash
# =====================================================================
# make_pas_index.sh - add PAS-related groups to index.ndx
#
# The peripheral anionic site (PAS) of human AChE (4EY6 numbering, mature
# protein) is formed by the aromatic residues lining the entrance of the
# gorge:
#
#     Tyr72, Asp74, Tyr124, Trp286, Tyr341
#
# The paper additionally identifies the 344-361 stretch as the peptide's main
# residence region on the AChE surface, and a new contact patch at Leu76/Trp77,
# so groups for those are created too:
#
#     PAS            Y72 D74 Y124 W286 Y341
#     PAS_aromatic   Y72 Y124 W286 Y341
#     Region344_361  residues 344-361      (main residence region)
#     Contact_L76_W77  L76 W77             (new contact patch)
#     CAS            W86 E202 S203 H447    (catalytic site, for reference)
#
# The residue numbers are taken from the RECEPTOR chain only, so a peptide
# residue with a clashing number is never picked up by mistake.
#
# Usage:
#   scripts/make_pas_index.sh -d work/alllhrc
#   scripts/make_pas_index.sh -d work/x --offset -4   # renumbered construct
# =====================================================================
set -euo pipefail

GMX="${GMX:-gmx}"
WORKDIR=""
STRUCT=""
NDX="index.ndx"
OFFSET=0

usage() {
    sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'
    cat <<EOF

Options:
  -d, --dir DIR      system directory (required)
  -s, --struct FILE  structure/tpr to use (default: md.tpr, em.tpr or solv_ions.gro)
  -n, --ndx FILE     index file to extend (default: $NDX)
      --offset N     add N to every PAS residue number (renumbered constructs)
  -h, --help         this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dir)    WORKDIR="$2"; shift 2 ;;
        -s|--struct) STRUCT="$2"; shift 2 ;;
        -n|--ndx)    NDX="$2"; shift 2 ;;
        --offset)    OFFSET="$2"; shift 2 ;;
        -h|--help)   usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -n "$WORKDIR" ]] || { echo "ERROR: -d/--dir is required" >&2; usage; exit 1; }
[[ -d "$WORKDIR" ]] || { echo "ERROR: no such directory: $WORKDIR" >&2; exit 1; }
cd "$WORKDIR"

if [[ -z "$STRUCT" ]]; then
    for c in md.tpr equil.tpr npt.tpr nvt.tpr em.tpr solv_ions.gro; do
        [[ -f "$c" ]] && { STRUCT="$c"; break; }
    done
fi
[[ -n "$STRUCT" && -f "$STRUCT" ]] || {
    echo "ERROR: no structure found; pass -s" >&2; exit 1; }
[[ -f "$NDX" ]] || { echo "ERROR: $NDX not found" >&2; exit 1; }

shift_res() { python3 -c "print(' '.join(str(int(r)+$OFFSET) for r in '$1'.split()))"; }

PAS_RES="$(shift_res "72 74 124 286 341")"
PAS_ARO="$(shift_res "72 124 286 341")"
REGION_LO="$(shift_res 344)"
REGION_HI="$(shift_res 361)"
PATCH="$(shift_res "76 77")"
CAS_RES="$(shift_res "86 202 203 447")"

echo "======================================================================"
echo " make_pas_index.sh"
echo "   structure : $STRUCT"
echo "   index     : $NDX"
echo "   offset    : $OFFSET"
echo "   PAS       : $PAS_RES"
echo "======================================================================"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

add_group() {
    local name="$1" sel="$2"
    local out="$TMP/${name}.ndx"
    if "$GMX" select -s "$STRUCT" -n "$NDX" -select "$sel" -on "$out" \
            >"$TMP/${name}.log" 2>&1; then
        local n
        n=$(grep -vc '^\[' "$out" 2>/dev/null || echo 0)
        # rename the group and append it
        {
            echo "[ $name ]"
            grep -v '^\[' "$out"
        } >> "$NDX"
        local natoms
        natoms=$(grep -v '^\[' "$out" | wc -w)
        printf "   %-18s %6s atoms\n" "$name" "$natoms"
    else
        # gmx select aborts on an empty selection, which almost always means
        # the residue numbers do not exist in THIS structure (renumbered
        # construct, or a receptor that is not full-length AChE)
        if grep -qi "never matches any atoms" "$TMP/${name}.log"; then
            echo "   WARNING: group '$name' is empty - those residue numbers"
            echo "            do not exist in the receptor. Check numbering"
            echo "            and use --offset if the construct is renumbered."
        else
            echo "   WARNING: could not build group '$name'"
            grep -iE "^(Error|Fatal|Inconsistency|Invalid)" "$TMP/${name}.log" \
                | head -3 | sed 's/^/       /'
        fi
    fi
}

# every selection is intersected with the Receptor group so peptide residues
# that happen to share a number are excluded
add_group "PAS"             "group \"Receptor\" and resid $PAS_RES"
add_group "PAS_aromatic"    "group \"Receptor\" and resid $PAS_ARO"
add_group "Region344_361"   "group \"Receptor\" and resid $REGION_LO to $REGION_HI"
add_group "Contact_L76_W77" "group \"Receptor\" and resid $PATCH"
add_group "CAS"             "group \"Receptor\" and resid $CAS_RES"

echo "----------------------------------------------------------------------"
echo " groups appended to $NDX"
echo " re-run scripts/analyze.sh to get mindist_peptide_PAS.xvg"
echo "======================================================================"
