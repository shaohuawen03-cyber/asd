#!/usr/bin/env bash
# ============================================================
# 建立用于分析的分组索引 (AChE / 肽 / 骨架 / 水)
# 请根据你的体系修改残基范围:
#   ACHERES : AChE 链的残基范围
#   PEPRES  : 肽链的残基范围 (pdb2gmx 后通常排在 AChE 之后)
# 用法: 在产物工作目录中运行
# ============================================================
set -euo pipefail

# AChE 残基范围 (请按实际体系修改, 4ey6 构建体约 1-537 或 1-548)
ACHERES="1-537"
# 肽残基范围 (AChE 之后的连续编号)
PEPRES="538-579"

gmx make_ndx -f neutral.gro -o index.ndx << EOF
ri ${ACHERES}
name 6 AChE
ri ${PEPRES}
name 7 Peptide
a CA C N O
name 8 Backbone
ri ${ACHERES} & a CA C N O
name 9 AChE_Backbone
ri ${PEPRES} & a CA C N O
name 10 Peptide_Backbone
keep 1
q
EOF
echo "已生成 index.ndx, 包含组: Protein(AChE+肽), AChE, Peptide, Backbone,"
echo "AChE_Backbone, Peptide_Backbone 等。请确认分组编号与残基范围正确。"
