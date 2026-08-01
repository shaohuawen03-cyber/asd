#!/usr/bin/env bash
# ============================================================
# 建立用于分析的分组索引 (AChE / 肽 / 骨架 / 水)
# 请根据你的体系修改残基范围:
#   ACHERES : AChE 链的残基范围
#   PEPRES  : 肽链的残基范围 (pdb2gmx 后通常排在 AChE 之后)
# 用法: 在产物工作目录中运行
# ============================================================
set -eu

if [ -n "${GMX:-}" ]; then
    gmx() { "${GMX}" "$@"; }
elif command -v gmx.exe >/dev/null 2>&1; then
    gmx() { gmx.exe "$@"; }
fi

# AChE 残基范围 (请按实际体系修改, 4ey6 构建体约 1-537)
ACHERES="${ACHERES:-1-537}"
# 肽残基范围 (AChE 之后的连续编号)
PEPRES="${PEPRES:-538-579}"

TPR_FILE="md.tpr"
if [ ! -f "${TPR_FILE}" ]; then
    TPR_FILE="neutral.gro"
fi

echo ">> 使用结构文件: ${TPR_FILE} 生成分析索引组 (ACHERES=${ACHERES}, PEPRES=${PEPRES}) ..."

gmx make_ndx -f "${TPR_FILE}" -o index.ndx << EOF
ri ${ACHERES}
name 19 AChE
ri ${PEPRES}
name 20 Peptide
19 & 4
name 21 AChE_Backbone
20 & 4
name 22 Peptide_Backbone
q
EOF

echo ">> 已成功生成 index.ndx, 包含新增分析组: AChE, Peptide, AChE_Backbone, Peptide_Backbone。"
