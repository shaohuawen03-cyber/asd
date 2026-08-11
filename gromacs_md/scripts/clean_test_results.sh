#!/usr/bin/env bash
# ============================================================
# Clean Up All Previous Test Results and Temporary Files (Bash / WSL)
#
# Usage:  ./clean_test_results.sh
# ============================================================
set -eu

DIRS=("../md_alllhrc" "../md_fllhttr" "../md_ylsllqr" "../md_ache" "../compare_figures")
FILES=("./test_four_systems.log")

echo "===================================================================="
echo " Cleaning up all previous test directories and log files..."
echo "===================================================================="

for D in "${DIRS[@]}"; do
    if [ -d "$D" ]; then
        rm -rf "$D"
        echo "  [DELETED DIRECTORY] $D"
    fi
done

for F in "${FILES[@]}"; do
    if [ -f "$F" ]; then
        rm -f "$F"
        echo "  [DELETED FILE] $F"
    fi
done

echo ""
echo "===================================================================="
echo " SUCCESS! Workspace is completely clean and ready for formal 100 ns simulations."
echo "===================================================================="
