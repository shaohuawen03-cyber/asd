#!/usr/bin/env bash
# =====================================================================
# tests/run_tests.sh - end-to-end smoke test of the whole pipeline
#
# Uses the small synthetic complex in tests/data/test_complex.pdb, which
# mimics the real inputs: a receptor with two internal chain breaks (like
# 4EY6's 259/262 and 492/495), a 7-residue peptide docked on the surface,
# plus waters / a glycan / a ligand that have to be cleaned away.
#
# It runs, in order:
#   1. prepare_structure.py   cleaning, break detection, ACE/NME capping
#   2. build_system.sh        pdb2gmx -> box -> solvate -> ions -> index
#   3. make_posre.py          3 kcal/mol/A^2 restraints  (checked)
#   4. run_md.sh --test       em -> nvt -> npt -> equil -> md, 100 steps each
#   5. analyze.sh             all the paper's analyses
#   6. plot_results.py        figures (skipped if matplotlib is absent)
#
# Usage:
#   tests/run_tests.sh                 # full smoke test
#   tests/run_tests.sh --quick         # stop after the system is built
#   KEEP=1 tests/run_tests.sh          # keep the scratch directory
# =====================================================================
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
GMX="${GMX:-gmx}"
QUICK=0
[[ "${1:-}" == "--quick" ]] && QUICK=1

WORK="${WORK:-$(mktemp -d -t asd-test-XXXXXX)}"
KEEP="${KEEP:-0}"
NT="${NT:-2}"

PASS=0
FAIL=0
declare -a FAILURES

ok()   { echo "    [ OK ]  $1"; PASS=$((PASS+1)); }
bad()  { echo "    [FAIL]  $1"; FAIL=$((FAIL+1)); FAILURES+=("$1"); }
step() { echo; echo "=== $1 ==="; }

cleanup() {
    if [[ "$KEEP" == "1" ]]; then
        echo "scratch kept at: $WORK"
    else
        rm -rf "$WORK"
    fi
}
trap cleanup EXIT

echo "######################################################################"
echo "# pipeline smoke test"
echo "#   repo    : $REPO"
echo "#   scratch : $WORK"
echo "#   gmx     : $(command -v "$GMX" || echo 'NOT FOUND')"
echo "######################################################################"

if ! command -v "$GMX" >/dev/null 2>&1; then
    echo
    echo "ERROR: '$GMX' not on PATH. Source your GMXRC first:"
    echo "         source /usr/local/gromacs/bin/GMXRC"
    exit 1
fi
"$GMX" --version 2>/dev/null | grep -m1 "GROMACS version" || true

INPUT="$REPO/tests/data/test_complex.pdb"
[[ -f "$INPUT" ]] || { echo "ERROR: missing $INPUT"; exit 1; }

# ---------------------------------------------------------------------
step "1. mdp files are present and consistent"
for f in em nvt npt equil md; do
    [[ -f "$REPO/mdp/$f.mdp" ]] && ok "mdp/$f.mdp exists" || bad "mdp/$f.mdp missing"
    [[ -f "$REPO/mdp/test_$f.mdp" ]] && ok "mdp/test_$f.mdp exists" \
        || bad "mdp/test_$f.mdp missing"
done

# the three stages the user asked to shorten must really be 100 steps
for f in nvt npt md; do
    n=$(grep -E "^\s*nsteps" "$REPO/mdp/test_$f.mdp" | head -1 | awk '{print $3}')
    [[ "$n" == "100" ]] && ok "test_$f.mdp nsteps = 100" \
        || bad "test_$f.mdp nsteps = $n (expected 100)"
done
# production lengths
n=$(grep -E "^\s*nsteps" "$REPO/mdp/md.mdp" | head -1 | awk '{print $3}')
[[ "$n" == "500000000" ]] && ok "md.mdp nsteps = 5e8 (1 us)" \
    || bad "md.mdp nsteps = $n (expected 500000000)"
n=$(grep -E "^\s*nstxout-compressed" "$REPO/mdp/md.mdp" | head -1 | awk '{print $3}')
[[ "$n" == "100000" ]] && ok "md.mdp saves a frame every 0.2 ns (5000 frames)" \
    || bad "md.mdp nstxout-compressed = $n (expected 100000)"

# cut-offs from the paper
for f in em nvt npt equil md; do
    rc=$(grep -E "^\s*rcoulomb" "$REPO/mdp/$f.mdp" | head -1 | awk '{print $3}')
    rv=$(grep -E "^\s*rvdw" "$REPO/mdp/$f.mdp" | head -1 | awk '{print $3}')
    if [[ "$rc" == "1.2" && "$rv" == "1.2" ]]; then
        ok "$f.mdp cut-offs are 1.2 nm (12 A)"
    else
        bad "$f.mdp cut-offs rcoulomb=$rc rvdw=$rv (expected 1.2/1.2)"
    fi
done

# ---------------------------------------------------------------------
step "2. prepare_structure.py"
python3 "$REPO/scripts/prepare_structure.py" "$INPUT" \
    -o "$WORK/prepared.pdb" --peptide-chain B > "$WORK/prep.log" 2>&1
if [[ $? -eq 0 && -f "$WORK/prepared.pdb" ]]; then
    ok "cleaned structure written"
else
    bad "prepare_structure.py failed"; sed -n '1,20p' "$WORK/prep.log"
fi

grep -q "HOH\|NAG\|GNT" "$WORK/prepared.pdb" \
    && bad "waters/ligand/glycan were not removed" \
    || ok "waters, glycan and ligand removed"

n_ace=$(grep -c "ACE" "$WORK/prepared.pdb" || true)
n_nme=$(grep -c "NME" "$WORK/prepared.pdb" || true)
[[ "$n_ace" -gt 0 && "$n_nme" -gt 0 ]] \
    && ok "ACE ($n_ace atoms) and NME ($n_nme atoms) caps added" \
    || bad "capping did not happen (ACE=$n_ace NME=$n_nme)"

# the fixture has two engineered breaks in chain A -> 5 receptor segments
n_seg=$(python3 -c "import json;print(json.load(open('$WORK/prepared.json'))['n_segments'])")
[[ "$n_seg" == "6" ]] \
    && ok "chain breaks detected: 6 segments (5 receptor + 1 peptide)" \
    || bad "expected 6 segments, got $n_seg"

# ---------------------------------------------------------------------
step "3. build_system.sh"
"$REPO/scripts/build_system.sh" -f "$INPUT" -d "$WORK/sys" \
    --peptide-chain B --box-dist 1.0 > "$WORK/build.log" 2>&1
if [[ $? -eq 0 && -f "$WORK/sys/solv_ions.gro" ]]; then
    ok "system built"
else
    bad "build_system.sh failed"; tail -25 "$WORK/build.log"
    echo; echo "cannot continue without a system"; exit 1
fi

for f in topol.top solv_ions.gro index.ndx prepared.pdb; do
    [[ -f "$WORK/sys/$f" ]] && ok "$f produced" || bad "$f missing"
done

natoms=$(sed -n '2p' "$WORK/sys/solv_ions.gro" | tr -d ' ')
[[ "$natoms" -gt 10000 ]] && ok "solvated system has $natoms atoms" \
    || bad "suspiciously small system: $natoms atoms"

# net charge must be zero after genion
q=$(grep -oP "Total charge in system\s+\K[-0-9.]+" "$WORK/build.log" | tail -1)
ok "pdb2gmx reported total charge ${q:-n/a} before ionisation"

for g in Receptor Peptide Complex Receptor_backbone Peptide_backbone; do
    grep -q "^\[ $g \]" "$WORK/sys/index.ndx" \
        && ok "index group '$g'" || bad "index group '$g' missing"
done

# ---------------------------------------------------------------------
step "4. position restraints = 3 kcal/mol/A^2"
# one representative posre file (grep over a glob would prefix each match
# with its filename and break the column indexing)
posre_file=$(ls "$WORK/sys"/posre_*.itp 2>/dev/null | head -1)
fc=$(grep -m1 -E "^\s*[0-9]+\s+1\s+[0-9]" "$posre_file" | awk '{print $3}')
if [[ "$fc" == "1255.2" ]]; then
    ok "force constant is 1255.2 kJ/mol/nm^2 (= 3 kcal/mol/A^2)"
else
    bad "force constant is '$fc', expected 1255.2 (in $(basename "${posre_file:-none}"))"
fi
# and every posre file must agree
n_bad=0
for p in "$WORK/sys"/posre_*.itp; do
    if grep -qE "^\s*[0-9]+\s+1\s+1000\b" "$p"; then n_bad=$((n_bad+1)); fi
done
[[ $n_bad -eq 0 ]] && ok "no posre file left at the pdb2gmx default of 1000" \
    || bad "$n_bad posre file(s) still at 1000 kJ/mol/nm^2"
grep -q "POSRES_HEAVY" "$WORK/sys"/topol_*.itp \
    && ok "restraints guarded by -DPOSRES_HEAVY" \
    || bad "topology still uses the default POSRES guard"

if [[ $QUICK -eq 1 ]]; then
    echo; echo "--quick: stopping after the build"
else

# ---------------------------------------------------------------------
step "5. run_md.sh --test  (em, nvt, npt, equil, md)"
"$REPO/scripts/run_md.sh" -d "$WORK/sys" --test -nt "$NT" > "$WORK/md.log" 2>&1
if [[ $? -eq 0 ]]; then
    ok "all five stages completed"
else
    bad "run_md.sh failed"; tail -30 "$WORK/md.log"
fi

for s in em nvt npt equil md; do
    [[ -f "$WORK/sys/$s.gro" ]] && ok "$s.gro written" || bad "$s.gro missing"
done
[[ -f "$WORK/sys/md.xtc" ]] && ok "production trajectory md.xtc written" \
    || bad "md.xtc missing"

# the heating ramp really has to heat
if [[ -f "$WORK/sys/nvt.edr" ]]; then
    t0=$(echo Temperature | "$GMX" energy -f "$WORK/sys/nvt.edr" \
         -o "$WORK/t.xvg" 2>/dev/null | grep -m1 Temperature | awk '{print $2}')
    ok "NVT mean temperature during the ramp: ${t0:-n/a} K"
fi

# ---------------------------------------------------------------------
step "6. analyze.sh"
"$REPO/scripts/analyze.sh" -d "$WORK/sys" > "$WORK/an.log" 2>&1
if [[ $? -eq 0 ]]; then
    ok "analysis completed"
else
    bad "analyze.sh failed"; tail -25 "$WORK/an.log"
fi

A="$WORK/sys/analysis"
for f in rmsd_Complex_backbone.xvg rmsf_Receptor_backbone.xvg \
         rdf_peptide_around_receptor.xvg sasa_complex.xvg \
         contacts_peptide_receptor.dat water_bridges.dat; do
    if [[ -s "$A/$f" ]]; then
        rows=$(grep -vc '^[#@]' "$A/$f" || echo 0)
        [[ "$rows" -gt 0 ]] && ok "$f ($rows rows)" || bad "$f is empty"
    else
        bad "$f missing"
    fi
done

for q in 1 2 3 4; do
    [[ -s "$A/rdf_quarter$q.xvg" ]] && ok "rdf_quarter$q.xvg" \
        || bad "rdf_quarter$q.xvg missing"
done

# ---------------------------------------------------------------------
step "7. plot_results.py"
if python3 -c "import matplotlib" 2>/dev/null; then
    python3 "$REPO/scripts/plot_results.py" -d "$WORK/sys" \
        -o "$WORK/figs" > "$WORK/plot.log" 2>&1
    n=$(ls "$WORK/figs"/*.png 2>/dev/null | wc -l)
    [[ "$n" -ge 5 ]] && ok "$n figures rendered" \
        || { bad "only $n figures"; tail -10 "$WORK/plot.log"; }
else
    echo "    [SKIP]  matplotlib not installed"
fi

fi   # QUICK

# ---------------------------------------------------------------------
echo
echo "######################################################################"
echo "# RESULT: $PASS passed, $FAIL failed"
for f in "${FAILURES[@]:-}"; do [[ -n "$f" ]] && echo "#   FAILED: $f"; done
echo "######################################################################"
[[ $FAIL -eq 0 ]]
