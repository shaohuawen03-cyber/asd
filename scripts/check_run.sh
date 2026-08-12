#!/usr/bin/env bash
# =====================================================================
# check_run.sh - is this run trustworthy?  And package what to share.
#
# Reads a finished (or in-progress) system directory and applies pass/fail
# criteria to the things that actually go wrong in an AChE-peptide MD:
#
#   1. did every stage finish, and how much of the production run is done
#   2. temperature / pressure / density / energy drift
#   3. did the peptide stay bound, or did it drift off into the solvent
#   4. is the receptor stable (backbone RMSD) and the peptide mobile
#   5. LINCS warnings, PME grid problems, other complaints in the logs
#   6. periodic-image distance - is the box still big enough
#
# It then writes a small tarball (a few hundred kB) with the numbers and
# logs needed for a second opinion, WITHOUT the huge trajectory files.
#
# Usage:
#   scripts/check_run.sh -d work/alllhrc
#   scripts/check_run.sh -d work/alllhrc --stage md --bundle
# =====================================================================
set -uo pipefail

GMX="${GMX:-gmx}"
WORKDIR=""
STAGE="md"
BUNDLE=0
OUT=""

usage() {
    sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'
    cat <<EOF

Options:
  -d, --dir DIR     system directory to check (required)
      --stage NAME  production stage to check (default: $STAGE)
      --bundle      also write a shareable tarball of the small files
  -o, --out FILE    tarball path (default: <dir>/../<name>_report.tar.gz)
  -h, --help        this help
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dir)   WORKDIR="$2"; shift 2 ;;
        --stage)    STAGE="$2"; shift 2 ;;
        --bundle)   BUNDLE=1; shift ;;
        -o|--out)   OUT="$2"; shift 2 ;;
        -h|--help)  usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

[[ -n "$WORKDIR" ]] || { echo "ERROR: -d/--dir is required" >&2; usage; exit 1; }
[[ -d "$WORKDIR" ]] || { echo "ERROR: no such directory: $WORKDIR" >&2; exit 1; }

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKDIR="$(cd "$WORKDIR" && pwd)"
NAME="$(basename "$WORKDIR")"
cd "$WORKDIR"

command -v "$GMX" >/dev/null 2>&1 || {
    echo "ERROR: '$GMX' not found. source your GMXRC or set GMX=..." >&2; exit 1; }

REPORT="$WORKDIR/check_report.txt"
: > "$REPORT"

PASS=0; WARN=0; FAIL=0
ok()   { echo "  [ OK ]   $1"   | tee -a "$REPORT"; PASS=$((PASS+1)); }
warn() { echo "  [WARN]   $1"   | tee -a "$REPORT"; WARN=$((WARN+1)); }
bad()  { echo "  [FAIL]   $1"   | tee -a "$REPORT"; FAIL=$((FAIL+1)); }
info() { echo "           $1"   | tee -a "$REPORT"; }
head2() { echo ""              | tee -a "$REPORT"
          echo "=== $1 ==="    | tee -a "$REPORT"; }

# numeric compare helpers (awk, so no bc dependency)
lt() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a+0 <  b+0)}'; }
gt() { awk -v a="$1" -v b="$2" 'BEGIN{exit !(a+0 >  b+0)}'; }

# mean of column 2 of an xvg
xvgmean() {
    [[ -s "$1" ]] || { echo ""; return; }
    grep -v '^[#@]' "$1" | awk '{s+=$2; n++} END{if(n)printf "%.3f", s/n}'
}
xvglast() {
    [[ -s "$1" ]] || { echo ""; return; }
    grep -v '^[#@]' "$1" | tail -1 | awk '{printf "%.3f", $2}'
}
xvgmax() {
    [[ -s "$1" ]] || { echo ""; return; }
    grep -v '^[#@]' "$1" | awk 'NR==1||$2>m{m=$2} END{printf "%.3f", m}'
}
# mean over the LAST n% of rows
xvgmean_tail() {
    [[ -s "$1" ]] || { echo ""; return; }
    local frac="$2"
    grep -v '^[#@]' "$1" | awk -v f="$frac" \
        '{v[NR]=$2} END{s=int(NR*(1-f)); if(s<1)s=1; t=0;c=0;
          for(i=s;i<=NR;i++){t+=v[i];c++} if(c)printf "%.3f", t/c}'
}
# energy term -> temporary xvg, echo the file path
energy_term() {
    local edr="$1" term="$2" out="$3"
    echo "$term" | "$GMX" energy -f "$edr" -o "$out" >/dev/null 2>&1
    [[ -s "$out" ]] && echo "$out" || echo ""
}

echo "######################################################################" | tee -a "$REPORT"
echo "# check_run.sh   -   $NAME" | tee -a "$REPORT"
echo "#   directory : $WORKDIR"   | tee -a "$REPORT"
echo "#   date      : $(date)"    | tee -a "$REPORT"
echo "#   gromacs   : $("$GMX" --version 2>/dev/null | awk '/GROMACS version/{print $3}')" | tee -a "$REPORT"
echo "######################################################################" | tee -a "$REPORT"

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# ---------------------------------------------------------------------
head2 "1. stages completed"
for s in em nvt npt equil "$STAGE"; do
    if [[ -f "$s.gro" ]]; then
        ok "$s finished ($s.gro present)"
    elif [[ -f "$s.cpt" ]]; then
        warn "$s started but has no final $s.gro (still running, or crashed)"
    else
        bad "$s did not run"
    fi
done

# how far did production actually get?
if [[ -f "$STAGE.log" ]]; then
    NSTEPS=$(grep -m1 -E "^\s*nsteps" "$STAGE.log" | awk '{print $3}')
    STEPDONE=$(grep -oE "^\s+Step\s+Time" -A1 "$STAGE.log" 2>/dev/null | tail -1 | awk '{print $1}')
    PERF=$(grep -m1 "Performance:" "$STAGE.log" | awk '{print $2}')
    [[ -n "${NSTEPS:-}" ]] && info "nsteps in the tpr : $NSTEPS"
    if [[ -n "${PERF:-}" ]]; then
        info "performance       : $PERF ns/day"
        if [[ -n "${NSTEPS:-}" ]]; then
            NS=$(awk -v n="$NSTEPS" 'BEGIN{printf "%.1f", n*0.002/1000}')
            DAYS=$(awk -v ns="$NS" -v p="$PERF" 'BEGIN{if(p>0)printf "%.1f", ns/p}')
            info "total length      : $NS ns  (~$DAYS days at this rate)"
        fi
    fi
fi

SHORT_RUN=0
if [[ -f "$STAGE.xtc" ]]; then
    NFR=$("$GMX" check -f "$STAGE.xtc" 2>&1 | tr '\r' '\n' | awk '/^Step/{print $2}' | tail -1)
    LAST=$("$GMX" check -f "$STAGE.xtc" 2>&1 | tr '\r' '\n' | awk '/Last frame/{print $NF}' | tail -1)
    if [[ -n "${NFR:-}" ]]; then
        info "trajectory        : $NFR frames, last time ${LAST:-?} ps"
        LASTNS=$(awk -v t="${LAST:-0}" 'BEGIN{printf "%.1f", t/1000}')
        if gt "$NFR" 4000; then
            ok "frame count $NFR - close to the paper's 5000"
        elif gt "$NFR" 400; then
            warn "only $NFR frames (${LASTNS} ns). Fine for a first look, but"
            info "         the paper averages over 5000 frames of a 1 us run."
        else
            SHORT_RUN=1
            warn "$NFR frames (${LASTNS} ns) - this looks like a TEST run."
            info "         Statistics from this are not publishable; the checks"
            info "         below are relaxed accordingly."
        fi
    fi
else
    bad "$STAGE.xtc missing - nothing to analyse"
fi

# ---------------------------------------------------------------------
head2 "2. thermodynamics"
if [[ -f "$STAGE.edr" ]]; then
    f=$(energy_term "$STAGE.edr" Temperature "$TMP/T.xvg")
    if [[ -n "$f" ]]; then
        T=$(xvgmean "$f"); TT=$(xvgmean_tail "$f" 0.5)
        info "mean temperature      : $T K   (last half: $TT K)"
        if gt "$T" 297 && lt "$T" 303; then
            ok "temperature is at the 300 K target"
        elif [[ $SHORT_RUN -eq 1 ]]; then
            warn "temperature $T K is off 300 K, but this is a short test run"
            info "         (a 100-step run never leaves the heating transient)"
        else
            bad "temperature $T K is off the 300 K target"
        fi
    fi
    f=$(energy_term "$STAGE.edr" Pressure "$TMP/P.xvg")
    if [[ -n "$f" ]]; then
        P=$(xvgmean "$f")
        info "mean pressure         : $P bar"
        # instantaneous pressure fluctuates by hundreds of bar; only the mean matters
        if gt "$P" -50 && lt "$P" 50; then ok "mean pressure is near 1 bar"
        elif [[ $SHORT_RUN -eq 1 ]]; then
            info "         mean pressure $P bar - meaningless over so few frames"
        else warn "mean pressure $P bar is far from 1 bar"; fi
    fi
    f=$(energy_term "$STAGE.edr" Density "$TMP/D.xvg")
    if [[ -n "$f" ]]; then
        D=$(xvgmean "$f"); DT=$(xvgmean_tail "$f" 0.5)
        info "mean density          : $D kg/m^3  (last half: $DT)"
        if gt "$D" 970 && lt "$D" 1050; then ok "density is physically sensible"
        else bad "density $D kg/m^3 is wrong - check solvation/ions"; fi
    fi
    f=$(energy_term "$STAGE.edr" Potential "$TMP/E.xvg")
    if [[ -n "$f" ]]; then
        E1=$(grep -v '^[#@]' "$f" | head -1 | awk '{print $2}')
        E2=$(grep -v '^[#@]' "$f" | tail -1 | awk '{print $2}')
        DR=$(awk -v a="$E1" -v b="$E2" 'BEGIN{if(a!=0)printf "%.3f", 100*(b-a)/(a<0?-a:a)}')
        info "potential energy      : $E1 -> $E2 kJ/mol  (${DR}% drift)"
        AB=$(awk -v d="$DR" 'BEGIN{printf "%.3f", (d<0?-d:d)}')
        if lt "$AB" 1.0; then ok "potential energy is stable (<1% drift)"
        elif [[ $SHORT_RUN -eq 1 ]]; then
            info "         ${DR}% drift - expected while the system is still heating"
        else warn "potential energy drifted ${DR}% - inspect $STAGE.log"; fi
    fi
else
    warn "$STAGE.edr missing - skipping the thermodynamic checks"
fi

# ---------------------------------------------------------------------
head2 "3. did the peptide stay bound?"
A="$WORKDIR/analysis"
MD_MIN="$A/mindist_peptide_receptor.xvg"
if [[ -s "$MD_MIN" ]]; then
    MEAN=$(xvgmean "$MD_MIN"); MX=$(xvgmax "$MD_MIN"); LASTV=$(xvglast "$MD_MIN")
    info "peptide-receptor min distance: mean $MEAN nm, max $MX nm, final $LASTV nm"
    if lt "$MEAN" 0.5; then
        ok "the peptide stayed in contact with AChE throughout"
    elif lt "$MEAN" 1.0; then
        warn "the peptide is loosely associated (mean $MEAN nm)"
    else
        bad "the peptide left the surface (mean $MEAN nm) - the complex dissociated"
    fi
    if gt "$MX" 1.5; then
        warn "at some point it was $MX nm away - it detached and may have re-bound"
    fi
else
    warn "no mindist_peptide_receptor.xvg - run scripts/analyze.sh first"
fi

PAS_MIN="$A/mindist_peptide_PAS.xvg"
if [[ -s "$PAS_MIN" ]]; then
    PM=$(xvgmean "$PAS_MIN")
    NEAR=$(grep -v '^[#@]' "$PAS_MIN" | awk '{n++; if($2<0.5)c++} END{if(n)printf "%.1f", 100*c/n}')
    info "peptide-PAS min distance     : mean $PM nm; within 5 A for ${NEAR}% of frames"
    if lt "$PM" 0.5; then ok "the peptide occupies the PAS"
    elif lt "$PM" 1.2; then
        warn "the peptide sits NEAR but not ON the PAS (mean $PM nm)"
        info "         The paper sees exactly this: Abeta migrates to the"
        info "         adjacent 344-361 region and only visits the PAS."
    else
        warn "the peptide is far from the PAS (mean $PM nm)"
    fi
else
    info "no PAS distance - run scripts/make_pas_index.sh then analyze.sh"
fi

# ---------------------------------------------------------------------
head2 "4. structural stability"
R_REC="$A/rmsd_Receptor_backbone.xvg"
if [[ -s "$R_REC" ]]; then
    M=$(xvgmean_tail "$R_REC" 0.5); MX=$(xvgmax "$R_REC")
    MA=$(awk -v v="$M" 'BEGIN{printf "%.2f", v*10}')
    MXA=$(awk -v v="$MX" 'BEGIN{printf "%.2f", v*10}')
    info "receptor backbone RMSD : mean(last half) $MA A, max $MXA A"
    if lt "$MA" 3.0; then ok "AChE is stable (paper: below 2 A)"
    elif lt "$MA" 5.0; then warn "AChE RMSD $MA A is higher than the paper's <2 A"
    else bad "AChE RMSD $MA A - the protein is falling apart"; fi
fi
R_PEP="$A/rmsd_Peptide_backbone.xvg"
if [[ -s "$R_PEP" ]]; then
    M=$(xvgmean_tail "$R_PEP" 0.5)
    MA=$(awk -v v="$M" 'BEGIN{printf "%.2f", v*10}')
    info "peptide backbone RMSD  : mean(last half) $MA A"
    if gt "$MA" 1.0; then
        ok "the peptide is mobile ($MA A) - expected for a short peptide"
    elif [[ $SHORT_RUN -eq 1 ]]; then
        info "         only $MA A - nothing moves in a 100-step run"
    else
        info "         the peptide barely moved ($MA A); unusual but not wrong"
    fi
fi
F_PEP="$A/rmsf_Peptide_backbone.xvg"
if [[ -s "$F_PEP" ]]; then
    MX=$(xvgmax "$F_PEP")
    MXA=$(awk -v v="$MX" 'BEGIN{printf "%.1f", v*10}')
    info "peptide max RMSF       : $MXA A  (paper reports 4-14 A for Abeta)"
fi

# ---------------------------------------------------------------------
head2 "5. warnings in the logs"
TOTW=0
for s in em nvt npt equil "$STAGE"; do
    [[ -f "$s.log" ]] || continue
    n=$(grep -ci "LINCS WARNING\|Water molecule.*can not be settled\|constraint.*deviation" "$s.log" 2>/dev/null | head -1)
    n=${n:-0}
    TOTW=$((TOTW+n))
    [[ "$n" -gt 0 ]] && warn "$s.log has $n constraint warnings"
done
[[ "$TOTW" -eq 0 ]] && ok "no LINCS/SETTLE warnings in any stage"

for s in "$STAGE"; do
    [[ -f "$s.log" ]] || continue
    if grep -qi "Fatal error" "$s.log"; then
        bad "$s.log contains a Fatal error"
        grep -A4 -i "Fatal error" "$s.log" | head -8 | sed 's/^/           /' | tee -a "$REPORT"
    fi
done

# pressure-coupling / box complaints
if [[ -f "$STAGE.log" ]] && grep -qi "Pressure scaling more than 1%" "$STAGE.log"; then
    warn "aggressive pressure scaling reported - box was still relaxing"
fi

# ---------------------------------------------------------------------
head2 "6. periodic images"
# molecules must be WHOLE for -pi to be meaningful; analyze.sh already made
# such a trajectory, otherwise build a throwaway one here
PITRAJ=""
if   [[ -f "${STAGE}_whole.xtc" ]]; then PITRAJ="${STAGE}_whole.xtc"
elif [[ -f "$STAGE.xtc" && -f "$STAGE.tpr" && -f index.ndx ]]; then
    if echo -e "Complex\nSystem" | "$GMX" trjconv -s "$STAGE.tpr" -f "$STAGE.xtc" \
            -n index.ndx -o "$TMP/whole.xtc" -pbc mol -center >/dev/null 2>&1; then
        PITRAJ="$TMP/whole.xtc"
    fi
fi
if [[ -n "$PITRAJ" && -f "$STAGE.tpr" && -f index.ndx ]]; then
    if echo "Complex" | "$GMX" mindist -s "$STAGE.tpr" -f "$PITRAJ" -n index.ndx \
            -pi -od "$TMP/pi.xvg" >"$TMP/pi.log" 2>&1; then
        MIND=$(grep -i "minimum distance" "$TMP/pi.log" | tail -1)
        MINPI=$(grep -v '^[#@]' "$TMP/pi.xvg" 2>/dev/null | awk 'NR==1||$2<m{m=$2} END{printf "%.3f", m}')
        if [[ -n "${MINPI:-}" ]]; then
            info "smallest periodic-image distance: $MINPI nm"
            if gt "$MINPI" 2.4; then ok "no periodic self-interaction (> 2x rvdw)"
            elif gt "$MINPI" 1.2; then warn "images came within $MINPI nm - marginal box size"
            else bad "images within $MINPI nm - the box was too small"; fi
        fi
    else
        info "(periodic-image check skipped)"
    fi
fi

# ---------------------------------------------------------------------
head2 "7. analysis outputs"
if [[ -d "$A" ]]; then
    n=$(ls "$A"/*.xvg "$A"/*.dat 2>/dev/null | wc -l)
    ok "$n analysis files in analysis/"
    for f in contacts_peptide_receptor.dat water_bridges.dat; do
        [[ -s "$A/$f" ]] && ok "$f present" || warn "$f missing"
    done
    if [[ -s "$A/contacts_peptide_receptor.dat" ]]; then
        echo "" | tee -a "$REPORT"
        info "top contacting peptide residues:"
        grep -v '^#' "$A/contacts_peptide_receptor.dat" | sort -k2 -gr | head -5 \
            | awk '{printf "             %-10s %8.1f contacts, occupancy %.0f%%\n",$1,$2,$5*100}' \
            | tee -a "$REPORT"
    fi
else
    warn "no analysis/ directory - run scripts/analyze.sh -d $WORKDIR"
fi

# ---------------------------------------------------------------------
echo ""                                                            | tee -a "$REPORT"
echo "######################################################################" | tee -a "$REPORT"
echo "# VERDICT: $PASS ok, $WARN warnings, $FAIL failures"          | tee -a "$REPORT"
if [[ $FAIL -gt 0 ]]; then
    echo "#   -> something is wrong; see the [FAIL] lines above."   | tee -a "$REPORT"
elif [[ $WARN -gt 0 ]]; then
    echo "#   -> usable, but read the [WARN] lines before trusting the numbers." | tee -a "$REPORT"
else
    echo "#   -> the run looks healthy."                            | tee -a "$REPORT"
fi
echo "######################################################################" | tee -a "$REPORT"
echo ""
echo "full report written to: $REPORT"

# ---------------------------------------------------------------------
if [[ $BUNDLE -eq 1 ]]; then
    OUT="${OUT:-$WORKDIR/../${NAME}_report.tar.gz}"
    BD="$TMP/${NAME}_report"
    mkdir -p "$BD/analysis" "$BD/logs" "$BD/mdout"

    cp "$REPORT" "$BD/" 2>/dev/null
    cp prepared.json "$BD/" 2>/dev/null
    [[ -d "$A" ]] && cp "$A"/*.xvg "$A"/*.dat "$A"/summary.txt "$BD/analysis/" 2>/dev/null

    # logs: full em/nvt/npt/equil, but only head+tail of the big md log
    for s in em nvt npt equil; do
        [[ -f "$s.log" ]] && cp "$s.log" "$BD/logs/" 2>/dev/null
    done
    if [[ -f "$STAGE.log" ]]; then
        { head -250 "$STAGE.log"; echo; echo "... [middle truncated] ..."; echo;
          tail -250 "$STAGE.log"; } > "$BD/logs/${STAGE}.log.trimmed"
    fi
    cp mdout_*.mdp "$BD/mdout/" 2>/dev/null
    # topology summary only - the full .top can be large and is regenerable
    if [[ -f topol.top ]]; then
        grep -v "^;" topol.top | grep -A100 "\[ molecules \]" > "$BD/molecules.txt" 2>/dev/null
    fi
    [[ -f index.ndx ]] && grep "^\[" index.ndx > "$BD/index_groups.txt" 2>/dev/null
    # first and last frame as PDB - tiny, and enough to eyeball the pose
    if [[ -f "$STAGE.tpr" && -f "$STAGE.xtc" && -f index.ndx ]]; then
        echo -e "Complex\nComplex" | "$GMX" trjconv -s "$STAGE.tpr" -f "$STAGE.xtc" \
            -n index.ndx -o "$BD/frame_first.pdb" -pbc mol -center -dump 0 >/dev/null 2>&1
        LASTT=$("$GMX" check -f "$STAGE.xtc" 2>&1 | tr '\r' '\n' | awk '/Last frame/{print $NF}' | tail -1)
        echo -e "Complex\nComplex" | "$GMX" trjconv -s "$STAGE.tpr" -f "$STAGE.xtc" \
            -n index.ndx -o "$BD/frame_last.pdb" -pbc mol -center -dump "${LASTT:-0}" >/dev/null 2>&1
    fi
    ( cd "$TMP" && tar czf "$OUT" "${NAME}_report" ) 2>/dev/null
    echo ""
    echo "======================================================================"
    echo " shareable bundle : $OUT"
    echo "                    $(du -h "$OUT" 2>/dev/null | cut -f1) - no trajectories inside"
    echo "======================================================================"
fi

[[ $FAIL -eq 0 ]]
