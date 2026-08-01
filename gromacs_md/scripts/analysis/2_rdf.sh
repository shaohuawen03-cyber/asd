#!/usr/bin/env bash
# ============================================================
# 图2: 肽(Aβ) 围绕 AChE 的径向分布函数 RDF
# 论文 3.1 节 (Aβ 围绕 AChE 在 3 A 内更易出现)
# gmx rdf 计算两质心之间的 RDF
# ============================================================
set -euo pipefail

# 整个产物轨迹的 RDF (图2A), 以质心计算 (mol_com)
gmx rdf -s md.tpr -f md.xtc -n index.ndx -o rdf_pep_ache.xvg \
       -ref AChE -sel Peptide -rdf mol_com -bin 0.02

# 将轨迹分四等份分别计算 RDF (图2B, 检验可靠性)
# 1000 ns 产物 => 每等份 250 ns
for q in 1 2 3 4; do
    t0=$(( (q-1) * 250 ))
    t1=$(( q * 250 ))
    gmx rdf -s md.tpr -f md.xtc -n index.ndx \
        -o rdf_pep_ache_q${q}.xvg -ref AChE -sel Peptide \
        -rdf mol_com -b ${t0} -e ${t1} -bin 0.02
done

echo "RDF 输出: rdf_pep_ache.xvg 及 rdf_pep_ache_q1..4.xvg"
