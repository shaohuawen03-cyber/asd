#!/usr/bin/env bash
# ============================================================
# 4-System 100 ns MD & Analysis Master Sequential Suite (User Protocol, Bash / WSL)
#
# Usage:  ./run_all_four_user_workflow.sh
# ============================================================
set -eu

SYSTEMS=("alllhrc" "fllhttr" "ylsllqr" "ache")
LOG_FILE="./run_all_four_user_workflow.log"

rm -f "$LOG_FILE"
unset TESTING || true
unset ONLY_MD || true

echo "===================================================================="
echo " [Step 0] Cleaning up previous test results to ensure 100% clean formal start..."
echo "===================================================================="
bash "./clean_test_results.sh"

echo ""
echo "===================================================================="
echo " Starting Sequential 100 ns MD & Analysis for 4 Systems (User Workflow)"
echo " Systems in queue: ${SYSTEMS[*]}"
echo " Log file: ${LOG_FILE}"
echo "===================================================================="

for SYS in "${SYSTEMS[@]}"; do
    echo ""
    echo "===================================================================="
    echo " >>> [STARTING FORMAL 100 NS SIMULATION & ANALYSIS] System: ${SYS}"
    echo "     Timestamp: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "===================================================================="

    echo "=== [100NS START] System: ${SYS} at $(date '+%Y-%m-%d %H:%M:%S') ===" >> "$LOG_FILE"
    bash "./run_split_md_workflow.sh" "${SYS}" 2>&1 | tee -a "$LOG_FILE"

    echo ""
    echo "===================================================================="
    echo " [READY FOR YOUR INSPECTION] System ${SYS} 100 ns MD and Figures Completed!"
    echo " -> Open directory: md_${SYS}/figures/ to inspect SVG/PNG/PDF plots & CSV tables!"
    echo " -> Moving directly to next system in queue without delay..."
    echo "===================================================================="
done

echo ""
echo "===================================================================="
echo " >>> [BONUS STEP] Generating AChE Monomer vs Complex Comparative Figures..."
echo "===================================================================="
for COMP_SYS in "alllhrc" "fllhttr" "ylsllqr"; do
    if [ -f "../md_ache/md_0_1.xtc" ] && [ -f "../md_${COMP_SYS}/md_0_1.xtc" ]; then
        echo "  -> Comparing ache Monomer vs ${COMP_SYS} Complex..."
        python3 "./analysis/plot_compare_systems.py" --protein "../md_ache" --complex "../md_${COMP_SYS}" --out "../compare_ache_vs_${COMP_SYS}" 2>&1 | tee -a "$LOG_FILE"
    fi
done

echo ""
echo "===================================================================="
echo " 4-SYSTEM 100 NS MD & ANALYSIS MASTER SUITE FINISHED!"
echo " Full execution log saved in: ${LOG_FILE}"
echo "===================================================================="
