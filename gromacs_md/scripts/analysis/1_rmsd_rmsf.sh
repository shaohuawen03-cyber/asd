#!/usr/bin/env bash
# ============================================================
# 图1: 复合物/AChE/肽 的骨架 C原子 RMSD 与 RMSF
# 论文 3.1 节
# ============================================================
set -euo pipefail

# 复合物骨架 RMSD (对齐到骨架, 论文图1A)
gmx rms -s md.tpr -f md.xtc -n index.ndx -o rmsd_complex_bb.xvg -fit rot+trans << EOF
Backbone
Backbone
EOF

# AChE 骨架 RMSD (论文图1C)
gmx rms -s md.tpr -f md.xtc -n index.ndx -o rmsd_ache_bb.xvg -fit rot+trans << EOF
AChE_Backbone
AChE_Backbone
EOF

# 肽 骨架 RMSD (论文图1E)
gmx rms -s md.tpr -f md.xtc -n index.ndx -o rmsd_pep_bb.xvg -fit rot+trans << EOF
Peptide_Backbone
Peptide_Backbone
EOF

# 复合物骨架 RMSF (论文图1B)
gmx rmsf -s md.tpr -f md.xtc -n index.ndx -o rmsf_complex_bb.xvg -res << EOF
Backbone
EOF

# AChE 骨架 RMSF (论文图1D)
gmx rmsf -s md.tpr -f md.xtc -n index.ndx -o rmsf_ache_bb.xvg -res << EOF
AChE_Backbone
EOF

# 肽 骨架 RMSF (论文图1F, 按残基)
gmx rmsf -s md.tpr -f md.xtc -n index.ndx -o rmsf_pep_bb.xvg -res << EOF
Peptide_Backbone
EOF

echo "RMSD/RMSF 输出: rmsd_*_bb.xvg, rmsf_*_bb.xvg"
echo "可用 xmgrace / Python matplotlib 作图, 或 gmx xvgconv 转换。"
