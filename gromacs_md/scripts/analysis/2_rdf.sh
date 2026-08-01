#!/usr/bin/env bash
# ============================================================
# 图2: 肽(Aβ) 围绕 AChE 的径向分布函数 RDF
# 论文 3.1 节 (Aβ 围绕 AChE 在 3 A 内更易出现)
# gmx rdf 计算两质心之间的 RDF
# ============================================================
set -eu

if [ -n "${GMX:-}" ]; then
    gmx() { "${GMX}" "$@"; }
elif command -v gmx.exe >/dev/null 2>&1; then
    gmx() { gmx.exe "$@"; }
fi

# 整个产物轨迹的 RDF (图2A), 以质心计算 (mol_com)
gmx rdf -s md.tpr -f md.xtc -n index.ndx -o rdf_pep_ache.xvg \
       -ref AChE -sel Peptide -selrpos mol_com -seltype mol_com -bin 0.02

# 将轨迹分四等份分别计算 RDF (图2B, 检验可靠性)
# 1000 ns 产物 => 每等份 250 ns (Gromacs -b/-e 单位为 ps -> 250000 ps)
# 测试模式 => 每等份 0.05 ps
TESTING="${TESTING:-0}"
if [ "${TESTING}" = "1" ]; then
    STEP_PS="0.05"
else
    STEP_PS="250000"
fi

for q in 1 2 3 4; do
    t0=$(awk "BEGIN {print ($q - 1) * $STEP_PS}")
    t1=$(awk "BEGIN {print $q * $STEP_PS}")
    gmx rdf -s md.tpr -f md.xtc -n index.ndx \
        -o "rdf_pep_ache_q${q}.xvg" -ref AChE -sel Peptide \
        -selrpos mol_com -seltype mol_com -b "${t0}" -e "${t1}" -bin 0.02 || echo ">> 跳过区间 Q${q} (测试轨迹跨度不足)"
done

echo "RDF 输出: rdf_pep_ache.xvg 及 rdf_pep_ache_q1..4.xvg"
