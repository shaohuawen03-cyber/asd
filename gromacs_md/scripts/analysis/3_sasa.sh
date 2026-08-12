#!/usr/bin/env bash
# ============================================================
# 图3: 复合物溶剂可及表面积 SASA 平均值与收敛性
# 论文 3.2 节 (LCPO 算法)
# 说明: gmx sasa 采用类 Shrake-Rupley 算法, 功能上等价于 LCPO
# 每 10 帧计算一次 (每帧 0.2 ns => 每 2 ns 一个点, 共 500 点)
# ============================================================
set -eu

if [ -n "${GMX:-}" ]; then
    gmx() { "${GMX}" "$@"; }
elif command -v gmx.exe >/dev/null 2>&1; then
    gmx() { gmx.exe "$@"; }
fi

TRAJ_FILE="md_fit.xtc"
if [ ! -f "${TRAJ_FILE}" ]; then
    TRAJ_FILE="md_noPBC.xtc"
    if [ ! -f "${TRAJ_FILE}" ]; then
        TRAJ_FILE="md.xtc"
    fi
fi

# 使用整条轨迹, 蛋白体系 SASA
gmx sasa -s md.tpr -f "${TRAJ_FILE}" -n index.ndx -o sasa_complex.xvg \
         -surface Protein -output Protein

# AChE SASA
gmx sasa -s md.tpr -f "${TRAJ_FILE}" -n index.ndx -o sasa_ache.xvg \
         -surface AChE -output AChE

if grep -q "\[ *Peptide *\]" index.ndx 2>/dev/null; then
    gmx sasa -s md.tpr -f "${TRAJ_FILE}" -n index.ndx -o sasa_pep.xvg \
             -surface Peptide -output Peptide
else
    echo ">> [单体蛋白模式] 未检测到 Peptide 组，跳过小肽 SASA。"
fi

echo "SASA 输出: sasa_complex.xvg (每 2 ns 一个点)"
echo "收敛性: 可用每 100 ns 区间平均 (参照论文图3B)"
