#!/usr/bin/env bash
# ============================================================
# AChE - A-beta Complex 100 ns Formal Production MD & Analysis Suite (Bash / WSL)
#
# Usage:  ./run_100ns_formal.sh alllhrc
#         ./run_100ns_formal.sh fllhttr
#         ./run_100ns_formal.sh ylsllqr
#         ./run_100ns_formal.sh ache
# ============================================================
set -eu

SYS="${1:?用法: ./run_100ns_formal.sh <系统前缀, 如 alllhrc>}"
WORK="../md_${SYS}"

echo "===================================================================="
echo " AChE(4ey6) - A-beta Complex 100 ns Formal Production MD & Analysis Pipeline"
echo " System: ${SYS}"
echo "===================================================================="

echo ""
echo ">>> [Step 0/3] Deleting any previous test results in ${WORK} to ensure a 100% clean start..."
if [ -d "${WORK}" ]; then
    rm -rf "${WORK}"
    echo "  [DELETED PREVIOUS DIRECTORY] ${WORK}"
fi

echo ""
echo ">>> [Step 1/3] Starting 100 ns Formal Production MD Simulation (mdp/100ns) ..."
bash "./run_all.sh" "${SYS}"

echo ""
echo ">>> [Step 2/3 & 3/3] Running Trajectory Analysis & Generating Publication SVG/PNG/PDF Figures ..."
bash "./run_analysis.sh" "${SYS}"

echo ""
echo "===================================================================="
echo " SUCCESS! System ${SYS} 100 ns Production MD simulation and analysis completed!"
echo " Publication Figures Directory: md_${SYS}/figures/"
echo "    - fig0_summary_all.{svg,png,pdf} (A-F 2x3 Combined Master Plot)"
echo "    - fig1_rmsd_rmsf to fig6_bridging_waters individual plots"
echo "    - summary_metrics.csv and wide.csv statistical summary tables"
echo "===================================================================="
