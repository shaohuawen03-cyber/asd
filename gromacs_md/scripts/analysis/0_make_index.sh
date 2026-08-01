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

# AChE / Peptide 选择表达式:
# 默认优先按 GROMACS 链ID (chain A = AChE, chain B = Peptide) 匹配，无论是 7 肽还是 42 肽皆能准确选中
# 若需手动指定残基范围，可在执行前传参: ACHE_SEL="ri 1-530" PEP_SEL="ri 531-537"
ACHE_SEL="${ACHE_SEL:-chain A}"
PEP_SEL="${PEP_SEL:-chain B}"

TPR_FILE="md.tpr"
if [ ! -f "${TPR_FILE}" ]; then
    TPR_FILE="neutral.gro"
fi

# 自动检测该体系 tpr 中原有默认组的最大编号(如 16 或 18 等)
LAST_IDX=$(echo "q" | gmx make_ndx -f "${TPR_FILE}" -o /dev/null 2>&1 | grep -E "^ *Group +[0-9]+" | tail -n 1 | awk '{print $2}' || true)
LAST_IDX="${LAST_IDX:-16}"

G1=$((LAST_IDX + 1))
G2=$((LAST_IDX + 2))
G3=$((LAST_IDX + 3))
G4=$((LAST_IDX + 4))

echo ">> 使用结构文件: ${TPR_FILE} 生成分析索引组 (ACHE_SEL='${ACHE_SEL}', PEP_SEL='${PEP_SEL}') ..."
echo ">> 原系统最大组号: ${LAST_IDX}, 自动分配新增分组编号: ${G1}(AChE), ${G2}(Peptide), ${G3}(AChE_Backbone), ${G4}(Peptide_Backbone) ..."

gmx make_ndx -f "${TPR_FILE}" -o index.ndx << EOF
${ACHE_SEL}
name ${G1} AChE
${PEP_SEL}
name ${G2} Peptide
${G1} & 4
name ${G3} AChE_Backbone
${G2} & 4
name ${G4} Peptide_Backbone
q
EOF

echo ">> 已成功生成 index.ndx, 包含新增分析组: AChE, Peptide, AChE_Backbone, Peptide_Backbone。"
