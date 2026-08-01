#!/usr/bin/env bash
# ============================================================
# 图4: 肽(Aβ) 二级结构倾向随时间的演变
# 论文 3.2 节 (helix, beta-turn, bend 倾向, 每 50 ns 统计)
# gmx do_dssp 基于 DSSP 判定二级结构, -ssdump 输出逐帧 SS 字符串
# ============================================================
set -euo pipefail

# 计算肽的二级结构随时间的分布
# ss_pep.sc : 每行一个帧, 形如 "<时间> <SS字符串>" (每残基一个 DSSP 码)
gmx do_dssp -s md.tpr -f md.xtc -n index.ndx \
            -ssdump ss_pep.sc -o ss_pep.xpm << EOF
Peptide
EOF

# 将 ssdump 输出分箱统计每 50 ns 的 helix / turn / bend 倾向
python3 dssp_bins.py ss_pep.sc ss_pep_bins.dat 50

echo "二级结构输出: ss_pep.xpm, ss_pep.sc, ss_pep_bins.dat"
echo "ss_pep_bins.dat 列: 时间窗口(ns)  helix倾向  turn倾向  bend倾向"
