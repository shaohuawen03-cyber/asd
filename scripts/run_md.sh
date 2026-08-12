#!/usr/bin/env bash
# =====================================================================
# run_md.sh - run the four-stage protocol of the reference paper
#
#   1. em     2000 steps minimisation,  heavy atoms restrained
#   2. nvt    1 ns heating 0 -> 300 K,  heavy atoms restrained
#   3. npt    1 ns density equilibration (NPT), heavy atoms restrained
#   4. equil  1 ns unrestrained equilibration
#   5. md     1000 ns production, frame every 0.2 ns -> 5000 frames
#
# --test swaps in mdp/test_*.mdp (100 steps each) so the whole chain can be
# verified in a couple of minutes.
#
# Usage:
#   scripts/run_md.sh -d work/alllhrc                # full production
#   scripts/run_md.sh -d work/alllhrc --test         # 100-step smoke test
#   scripts/run_md.sh -d work/alllhrc --stage md     # one stage only
#   scripts/run_md.sh -d work/alllhrc -nt 16 --gpu   # resources
# =====================================================================
set -euo pipefail

GMX="${GMX:-gmx}"
WORKDIR=""
TEST=0
STAGE="all"
NT=""
NTMPI=""
GPU=0
MAXWARN=2
MDRUN_EXTRA=""

usage() {
    sed -n '2,18p' "$0" | sed 's/^# \{0,1\}//'
    cat <<EOF

Options:
  -d, --dir DIR        system directory made by build_system.sh (required)
      --test           use the short mdp/test_*.mdp files (100 steps)
      --stage NAME     em | nvt | npt | equil | md | all   (default: all)
  -nt, --nt N          total threads for mdrun
      --ntmpi N        thread-MPI ranks
      --gpu            add -nb gpu -pme gpu -bonded gpu
      --maxwarn N      grompp -maxwarn (default: $MAXWARN)
      --mdrun-extra S  extra flags passed verbatim to mdrun
  -h, --help           this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dir)        WORKDIR="$2"; shift 2 ;;
        --test)          TEST=1; shift ;;
        --stage)         STAGE="$2"; shift 2 ;;
        -nt|--nt)        NT="$2"; shift 2 ;;
        --ntmpi)         NTMPI="$2"; shift 2 ;;
        --gpu)           GPU=1; shift ;;
        --maxwarn)       MAXWARN="$2"; shift 2 ;;
        --mdrun-extra)   MDRUN_EXTRA="$2"; shift 2 ;;
        -h|--help)       usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -n "$WORKDIR" ]] || { echo "ERROR: -d/--dir is required" >&2; usage; exit 1; }
[[ -d "$WORKDIR" ]] || { echo "ERROR: no such directory: $WORKDIR" >&2; exit 1; }

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MDPDIR="$REPO/mdp"
WORKDIR="$(cd "$WORKDIR" && pwd)"
cd "$WORKDIR"

command -v "$GMX" >/dev/null 2>&1 || {
    echo "ERROR: '$GMX' not found. source your GMXRC or set GMX=..." >&2; exit 1; }

[[ -f topol.top      ]] || { echo "ERROR: topol.top missing - run build_system.sh" >&2; exit 1; }
[[ -f solv_ions.gro  ]] || { echo "ERROR: solv_ions.gro missing - run build_system.sh" >&2; exit 1; }

PREFIX=""
LABEL="PRODUCTION"
if [[ $TEST -eq 1 ]]; then
    PREFIX="test_"
    LABEL="TEST (short)"
fi

# ------------------------------------------------------------ mdrun flags
MDRUN_ARGS=()
[[ -n "$NT"    ]] && MDRUN_ARGS+=(-nt "$NT")
[[ -n "$NTMPI" ]] && MDRUN_ARGS+=(-ntmpi "$NTMPI")
if [[ $GPU -eq 1 ]]; then
    MDRUN_ARGS+=(-nb gpu -pme gpu -bonded gpu)
fi
# shellcheck disable=SC2206
[[ -n "$MDRUN_EXTRA" ]] && MDRUN_ARGS+=($MDRUN_EXTRA)

NDX_ARGS=()
[[ -f index.ndx ]] && NDX_ARGS+=(-n index.ndx)

echo "======================================================================"
echo " run_md.sh   [$LABEL]"
echo "   workdir : $WORKDIR"
echo "   mdp dir : $MDPDIR  (prefix '${PREFIX}')"
echo "   stage   : $STAGE"
echo "======================================================================"

# ------------------------------------------------------------ stage runner
# $1 stage name   $2 input coordinates   $3 restraint-reference coordinates
run_stage() {
    local name="$1" struct="$2" ref="$3"
    local mdp="$MDPDIR/${PREFIX}${name}.mdp"

    [[ -f "$mdp" ]] || { echo "ERROR: missing $mdp" >&2; exit 1; }
    [[ -f "$struct" ]] || { echo "ERROR: missing input coordinates $struct" >&2; exit 1; }

    echo
    echo "----------------------------------------------------------------------"
    echo ">>> STAGE $name   ($(basename "$mdp"))"
    echo "----------------------------------------------------------------------"

    local grompp_args=(-f "$mdp" -c "$struct" -p topol.top -o "${name}.tpr"
                       -po "mdout_${name}.mdp" -maxwarn "$MAXWARN")
    # -r is required whenever position restraints are active
    if grep -qE '^\s*define\s*=.*POSRES' "$mdp"; then
        grompp_args+=(-r "$ref")
    fi
    # velocities are inherited from the previous stage via the checkpoint
    if [[ "$name" != "em" && "$name" != "nvt" ]]; then
        local prev_cpt
        case "$name" in
            npt)   prev_cpt="nvt.cpt" ;;
            equil) prev_cpt="npt.cpt" ;;
            md)    prev_cpt="equil.cpt" ;;
        esac
        [[ -f "$prev_cpt" ]] && grompp_args+=(-t "$prev_cpt")
    fi
    grompp_args+=("${NDX_ARGS[@]}")

    "$GMX" grompp "${grompp_args[@]}"

    local mdrun_args=(-deffnm "$name" -v "${MDRUN_ARGS[@]}")
    # resume automatically if a checkpoint for THIS stage already exists
    if [[ -f "${name}.cpt" ]]; then
        echo ">>> found ${name}.cpt - resuming"
        mdrun_args+=(-cpi "${name}.cpt")
    fi
    "$GMX" mdrun "${mdrun_args[@]}"

    echo ">>> stage $name finished -> ${name}.gro"
}

case "$STAGE" in
    em)    run_stage em    solv_ions.gro solv_ions.gro ;;
    nvt)   run_stage nvt   em.gro        em.gro ;;
    npt)   run_stage npt   nvt.gro       nvt.gro ;;
    equil) run_stage equil npt.gro       npt.gro ;;
    md)    run_stage md    equil.gro     equil.gro ;;
    all)
        run_stage em    solv_ions.gro solv_ions.gro
        run_stage nvt   em.gro        em.gro
        run_stage npt   nvt.gro       nvt.gro
        run_stage equil npt.gro       npt.gro
        run_stage md    equil.gro     equil.gro
        ;;
    *) echo "ERROR: unknown stage '$STAGE'" >&2; exit 1 ;;
esac

echo
echo "======================================================================"
echo " DONE [$LABEL]"
if [[ -f md.xtc ]]; then
    echo "   trajectory : $WORKDIR/md.xtc"
    echo "   run input  : $WORKDIR/md.tpr"
    echo
    echo " Next:  scripts/analyze.sh -d $WORKDIR"
fi
echo "======================================================================"
