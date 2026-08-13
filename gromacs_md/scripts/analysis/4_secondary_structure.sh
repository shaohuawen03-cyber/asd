#!/usr/bin/env bash
# ============================================================
# Fig 4: COMPLEX (AChE+peptide) DSSP is the main product.
# Peptide-only DSSP is kept as a small supplement (7-mer is mostly coil).
# ============================================================
set -eu

if [ -n "${GMX:-}" ]; then
    gmx() { "${GMX}" "$@"; }
elif command -v gmx.exe >/dev/null 2>&1; then
    gmx() { gmx.exe "$@"; }
fi

for ana_dir in "/mnt/f/anaconda" "/mnt/f/Anaconda" "/f/anaconda" "/f/Anaconda" "/mnt/c/anaconda3" "/c/anaconda3" "$HOME/anaconda3" "$HOME/Anaconda3"; do
    if [ -d "$ana_dir" ]; then
        export PATH="$ana_dir:$ana_dir/Scripts:$ana_dir/Library/bin:$PATH"
        break
    fi
done

PY_FOUND=""
if [ -n "${PYTHON:-}" ]; then
    PY_FOUND="${PYTHON}"
else
    for cand in \
        "/mnt/f/anaconda/python.exe" \
        "/mnt/f/Anaconda/python.exe" \
        "/mnt/c/anaconda3/python.exe" \
        "/f/anaconda/python.exe" \
        "/f/Anaconda/python.exe" \
        "/c/anaconda3/python.exe" \
        "python.exe" \
        "python" \
        "python3"; do
        if command -v "$cand" >/dev/null 2>&1 || [ -x "$cand" ]; then
            if "$cand" -c "import numpy" >/dev/null 2>&1; then
                PY_FOUND="$cand"
                break
            fi
        fi
    done
fi
PY="${PY_FOUND:-python}"

TPR_FILE="md_0_1.tpr"
if [ ! -f "${TPR_FILE}" ]; then
    TPR_FILE="md.tpr"
fi
TRAJ_FILE="md_fit.xtc"
if [ ! -f "${TRAJ_FILE}" ]; then
    TRAJ_FILE="md_0_1.xtc"
    if [ ! -f "${TRAJ_FILE}" ]; then
        TRAJ_FILE="md.xtc"
    fi
fi

run_gmx_dssp() {
    local sel="$1"
    local odat="$2"
    local onum="$3"
    echo ">> [gmx dssp] -sel ${sel}  -> ${odat}"
    gmx dssp -s "${TPR_FILE}" -f "${TRAJ_FILE}" -n index.ndx \
        -sel "${sel}" -o "${odat}" -num "${onum}" \
        -hmode dssp -clear
}

set +e
if gmx dssp -h 2>&1 | grep -q -- "-sel "; then
    echo ">> PRIMARY: complex / Protein DSSP (AChE + peptide, ~530 residues)"
    for SEL in 'group "Protein" or group "Peptide"' 'group Protein or group Peptide' 'group Protein' 'resid 1 to 537'; do
        if run_gmx_dssp "${SEL}" ss_complex.dat ss_complex_num.xvg; then
            NROW=$(grep -v '^[#@;]' ss_complex_num.xvg 2>/dev/null | grep -c '[0-9]' || true)
            SLEN=$(grep -v '^[#@;]' ss_complex.dat 2>/dev/null | head -n 1 | tr -d '[:space:]' | wc -c || true)
            echo ">> complex dssp frames=${NROW}  SS_len=${SLEN}"
            if [ "${NROW}" -ge 10 ] && [ "${SLEN}" -ge 50 ]; then
                break
            fi
        fi
    done
    if grep -q "\[ *Peptide *\]" index.ndx 2>/dev/null; then
        echo ">> SUPPLEMENT: peptide-only DSSP"
        run_gmx_dssp 'group "Peptide"' ss_pep.dat ss_pep_num.xvg || \
            run_gmx_dssp "Peptide" ss_pep.dat ss_pep_num.xvg || true
    fi
fi
set -e

TESTING="${TESTING:-0}"
if [ "${TESTING}" = "1" ]; then
    WIN_NS=0.05
else
    WIN_NS=1
fi

to_py_path() {
    local p="$1"
    if command -v wslpath >/dev/null 2>&1 && [[ "${PY}" == *".exe"* || "${PY}" == *"/mnt/"* || "${PY}" == *":"* ]]; then
        wslpath -w "$p" 2>/dev/null || echo "$p"
    else
        echo "$p"
    fi
}

SS_SCRIPT="$(dirname "$0")/compute_peptide_ss.py"
echo ">> [DSSP rebuild] complex FIRST, then peptide; fix 500 ns fake time axis"
"${PY}" "$(to_py_path "${SS_SCRIPT}")" --dir . --nres 7 --window-ns "${WIN_NS}" --target both || true
echo "outputs: ss_complex_frac.xvg / ss_complex_summary.txt / ss_pep_frac.xvg"
