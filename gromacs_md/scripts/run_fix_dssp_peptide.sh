#!/usr/bin/env bash
# 只重跑有问题的部分: 肽 DSSP + 肽 RMSD 三段相位 + 重新出图
# 用法: ./run_fix_dssp_peptide.sh alllhrc
set -eu

SYS="${1:?用法: ./run_fix_dssp_peptide.sh <前缀, 如 alllhrc>}"
WORK="../md_${SYS}"

if [ ! -d "${WORK}" ]; then
    echo "!!! 错误: 找不到 ${WORK}"
    exit 1
fi

if [ -n "${PYTHON:-}" ]; then
    PY="${PYTHON}"
elif command -v python3 >/dev/null 2>&1; then
    PY="python3"
else
    PY="python"
fi

echo "========== 只修复 DSSP + 肽三段 RMSD: ${SYS} =========="
cd "${WORK}"

if [ ! -f index.ndx ]; then
    echo "[0] 补建 index.ndx ..."
    bash "../scripts/analysis/0_make_index.sh"
else
    echo "[0] 已有 index.ndx, 跳过 make_ndx / trjconv"
fi

echo "[1/3] 肽 DSSP 逐帧重建 ..."
bash "../scripts/analysis/4_secondary_structure.sh" || true

echo "[2/4] Peptide RMSD jump diagnosis ..."
"${PY}" "../scripts/analysis/diagnose_peptide_rmsd.py" -d . || true

echo "[3/4] Optional phase table ..."
"${PY}" "../scripts/analysis/analyze_peptide_phases.py" -d . || true

echo "[4/4] Replot figures ..."
"${PY}" "../scripts/analysis/plot_all.py" --dir . --out ./figures

echo "========== 完成: ${WORK}/figures/ =========="
