#!/usr/bin/env bash
# =====================================================================
# run_all.sh - build, simulate and analyse every complex in input/
#
# Runs the whole protocol for each *_complex.pdb found, one after another,
# then produces a side-by-side comparison of the peptides.
#
# Usage:
#   scripts/run_all.sh --test                 # smoke test, all systems
#   scripts/run_all.sh                        # full 1 us production
#   scripts/run_all.sh -i input -w work -nt 16 --gpu
# =====================================================================
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INDIR="$REPO/input"
WORKROOT="$REPO/work"
TEST=0
NT=""
GPU=""
PEPTIDE_CHAIN=""
BUILD_ONLY=0
SKIP_ANALYSIS=0
EXTRA_BUILD=""

usage() {
    sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'
    cat <<EOF

Options:
  -i, --input DIR      directory with *_complex.pdb (default: $INDIR)
  -w, --work DIR       root for the per-system directories (default: $WORKROOT)
      --test           100-step smoke test instead of the production run
      --peptide-chain  chain ID of the peptide (passed to build_system.sh)
  -nt, --nt N          threads for mdrun
      --gpu            enable GPU offload
      --build-only     stop after building the systems
      --skip-analysis  simulate but do not analyse
      --build-extra S  extra flags for build_system.sh
  -h, --help           this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)       INDIR="$2"; shift 2 ;;
        -w|--work)        WORKROOT="$2"; shift 2 ;;
        --test)           TEST=1; shift ;;
        --peptide-chain)  PEPTIDE_CHAIN="$2"; shift 2 ;;
        -nt|--nt)         NT="$2"; shift 2 ;;
        --gpu)            GPU="--gpu"; shift ;;
        --build-only)     BUILD_ONLY=1; shift ;;
        --skip-analysis)  SKIP_ANALYSIS=1; shift ;;
        --build-extra)    EXTRA_BUILD="$2"; shift 2 ;;
        -h|--help)        usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -d "$INDIR" ]] || { echo "ERROR: no such input directory: $INDIR" >&2; exit 1; }

shopt -s nullglob
PDBS=("$INDIR"/*.pdb)
shopt -u nullglob
[[ ${#PDBS[@]} -gt 0 ]] || { echo "ERROR: no .pdb files in $INDIR" >&2; exit 1; }

mkdir -p "$WORKROOT"

TEST_FLAG=""
[[ $TEST -eq 1 ]] && TEST_FLAG="--test"
NT_FLAG=()
[[ -n "$NT" ]] && NT_FLAG=(-nt "$NT")

echo "######################################################################"
echo "# run_all.sh"
echo "#   inputs : ${#PDBS[@]} structures from $INDIR"
echo "#   work   : $WORKROOT"
echo "#   mode   : $([[ $TEST -eq 1 ]] && echo 'TEST (100 steps)' || echo 'PRODUCTION (1 us)')"
echo "######################################################################"

DONE=()
FAILED=()

for pdb in "${PDBS[@]}"; do
    name="$(basename "$pdb" .pdb)"
    name="${name%_complex}"
    wd="$WORKROOT/$name"

    echo
    echo "######################################################################"
    echo "# SYSTEM: $name"
    echo "######################################################################"

    build_args=(-f "$pdb" -d "$wd")
    [[ -n "$PEPTIDE_CHAIN" ]] && build_args+=(--peptide-chain "$PEPTIDE_CHAIN")
    # shellcheck disable=SC2206
    [[ -n "$EXTRA_BUILD" ]] && build_args+=($EXTRA_BUILD)

    if [[ -f "$wd/solv_ions.gro" ]]; then
        echo ">>> already built, skipping build_system.sh"
    elif ! "$REPO/scripts/build_system.sh" "${build_args[@]}"; then
        echo "!!! BUILD FAILED for $name" >&2
        FAILED+=("$name (build)")
        continue
    fi

    # PAS groups are only meaningful for a full-length AChE receptor; a failure
    # here is not fatal for the rest of the pipeline
    "$REPO/scripts/make_pas_index.sh" -d "$wd" || true

    [[ $BUILD_ONLY -eq 1 ]] && { DONE+=("$name"); continue; }

    if ! "$REPO/scripts/run_md.sh" -d "$wd" $TEST_FLAG "${NT_FLAG[@]}" $GPU; then
        echo "!!! MD FAILED for $name" >&2
        FAILED+=("$name (md)")
        continue
    fi

    if [[ $SKIP_ANALYSIS -eq 0 ]]; then
        if ! "$REPO/scripts/analyze.sh" -d "$wd"; then
            echo "!!! ANALYSIS FAILED for $name" >&2
            FAILED+=("$name (analysis)")
            continue
        fi
    fi

    DONE+=("$name")
done

# ---------------------------------------------------------------- comparison
if [[ ${#DONE[@]} -gt 0 && $BUILD_ONLY -eq 0 && $SKIP_ANALYSIS -eq 0 ]]; then
    echo
    echo "######################################################################"
    echo "# comparison figures"
    echo "######################################################################"
    dirs=()
    for n in "${DONE[@]}"; do dirs+=("$WORKROOT/$n"); done
    python3 "$REPO/scripts/plot_results.py" --compare "${dirs[@]}" \
        -o "$WORKROOT/comparison" || \
        echo "    (plotting skipped - matplotlib not installed?)"
fi

echo
echo "######################################################################"
echo "# SUMMARY"
for n in "${DONE[@]}";   do echo "#   OK      $n"; done
for n in "${FAILED[@]}"; do echo "#   FAILED  $n"; done
echo "######################################################################"
[[ ${#FAILED[@]} -eq 0 ]]
