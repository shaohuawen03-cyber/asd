#!/usr/bin/env bash
# Radius of gyration (complex / AChE). Inserted before H-bond analysis.
set -eu

if [ -n "${GMX:-}" ]; then
    gmx() { "${GMX}" "$@"; }
elif command -v gmx.exe >/dev/null 2>&1; then
    gmx() { gmx.exe "$@"; }
fi

TPR_FILE="md.tpr"
if [ ! -f "${TPR_FILE}" ]; then
    TPR_FILE="md_0_1.tpr"
fi
TRAJ_FILE="md_fit.xtc"
if [ ! -f "${TRAJ_FILE}" ]; then
    TRAJ_FILE="md.xtc"
    if [ ! -f "${TRAJ_FILE}" ]; then
        TRAJ_FILE="md_0_1.xtc"
    fi
fi

if grep -q "\[ *Protein *\]" index.ndx 2>/dev/null; then
    gmx gyrate -s "${TPR_FILE}" -f "${TRAJ_FILE}" -n index.ndx -o gyrate_complex.xvg << EOF
Protein
EOF
else
    gmx gyrate -s "${TPR_FILE}" -f "${TRAJ_FILE}" -n index.ndx -o gyrate_complex.xvg << EOF
Backbone
EOF
fi

if grep -q "\[ *AChE *\]" index.ndx 2>/dev/null; then
    gmx gyrate -s "${TPR_FILE}" -f "${TRAJ_FILE}" -n index.ndx -o gyrate_ache.xvg << EOF
AChE
EOF
fi

echo "Rg output: gyrate_complex.xvg (and gyrate_ache.xvg if AChE group exists)"
