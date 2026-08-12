#!/usr/bin/env bash
# ============================================================
# 图1: 复合物/AChE/肽 的骨架 C原子 RMSD 与 RMSF
# 论文 3.1 节
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

# 复合物骨架 RMSD (对齐到骨架, 论文图1A)
gmx rms -s md.tpr -f "${TRAJ_FILE}" -n index.ndx -o rmsd_complex_bb.xvg -fit rot+trans << EOF
Backbone
Backbone
EOF

# AChE 骨架 RMSD (论文图1C)
gmx rms -s md.tpr -f "${TRAJ_FILE}" -n index.ndx -o rmsd_ache_bb.xvg -fit rot+trans << EOF
AChE_Backbone
AChE_Backbone
EOF

# 复合物骨架 RMSF (论文图1B)
gmx rmsf -s md.tpr -f "${TRAJ_FILE}" -n index.ndx -o rmsf_complex_bb.xvg -res << EOF
Backbone
EOF

# AChE 骨架 RMSF (论文图1D)
gmx rmsf -s md.tpr -f "${TRAJ_FILE}" -n index.ndx -o rmsf_ache_bb.xvg -res << EOF
AChE_Backbone
EOF

# 肽 骨架 RMSD/RMSF (复合物体系执行, 单体体系跳过)
if grep -q "\[ *Peptide_Backbone *\]" index.ndx 2>/dev/null; then
    gmx rms -s md.tpr -f "${TRAJ_FILE}" -n index.ndx -o rmsd_pep_bb.xvg -fit rot+trans << EOF
Peptide_Backbone
Peptide_Backbone
EOF
    if grep -q "\[ *Peptide_Calpha *\]" index.ndx 2>/dev/null; then
        gmx rmsf -s md.tpr -f "${TRAJ_FILE}" -n index.ndx -o rmsf_pep_bb.xvg -res << EOF
Peptide_Calpha
EOF
    else
        gmx rmsf -s md.tpr -f "${TRAJ_FILE}" -n index.ndx -o rmsf_pep_bb.xvg -res << EOF
Peptide_Backbone
EOF
    fi
else
    echo ">> [单体蛋白模式] 未检测到 Peptide_Backbone，跳过肽相关 RMSD / RMSF。"
fi

echo "RMSD/RMSF 输出: rmsd_*_bb.xvg, rmsf_*_bb.xvg"
echo "可用 xmgrace / Python matplotlib 作图, 或 gmx xvgconv 转换。"
