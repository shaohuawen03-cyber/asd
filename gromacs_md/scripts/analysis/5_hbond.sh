#!/usr/bin/env bash
# ============================================================
# 论文 3.3 节: 氢键分析
#   - AChE 与肽之间的氢键 (总 164 个)
#   - 肽内部氢键 (总 57383 个)
# gmx hbond 按几何判据 (距离截断 3.0 A, 角度 30 deg)
# ============================================================
set -eu

if [ -n "${GMX:-}" ]; then
    gmx() { "${GMX}" "$@"; }
elif command -v gmx.exe >/dev/null 2>&1; then
    gmx() { gmx.exe "$@"; }
fi

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

if ! grep -q "\[ *Peptide *\]" index.ndx 2>/dev/null; then
    echo ">> [单体蛋白模式] 未检测到 Peptide 组，跳过结合氢键计算，仅计算 AChE 内部氢键。"
    gmx hbond -s "${TPR_FILE}" -f "${TRAJ_FILE}" -n index.ndx -num hbond_ache_intra.xvg << EOF
AChE
AChE
EOF
    exit 0
fi

# AChE 与肽之间氢键
gmx hbond -s "${TPR_FILE}" -f "${TRAJ_FILE}" -n index.ndx -num hbond_ache_pep.xvg << EOF
AChE
Peptide
EOF

# 肽内部氢键
gmx hbond -s "${TPR_FILE}" -f "${TRAJ_FILE}" -n index.ndx -num hbond_pep_intra.xvg << EOF
Peptide
Peptide
EOF

# AChE 内部氢键 (可选)
gmx hbond -s "${TPR_FILE}" -f "${TRAJ_FILE}" -n index.ndx -num hbond_ache_intra.xvg << EOF
AChE
AChE
EOF

# 逐残基氢键寿命/占有率可加 -life, -dist, -ang 等选项
echo "氢键输出: hbond_ache_pep.xvg, hbond_pep_intra.xvg, hbond_ache_intra.xvg"
echo "xvg 第一列为时间, 第二列为该时刻氢键数; 统计其总和/平均可复现论文 3.3。"
