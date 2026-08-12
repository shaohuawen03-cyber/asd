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

if ! grep -q "\[ *Peptide *\]" index.ndx 2>/dev/null; then
    echo ">> [单体蛋白模式] 未检测到 Peptide 组，跳过肽-AChE 径向分布函数 RDF 计算。"
    exit 0
fi

TRAJ_FILE="md_fit.xtc"
if [ ! -f "${TRAJ_FILE}" ]; then
    TRAJ_FILE="md_noPBC.xtc"
    if [ ! -f "${TRAJ_FILE}" ]; then
        TRAJ_FILE="md.xtc"
    fi
fi

# 整个产物轨迹的 RDF (图2A), 以质心计算 (mol_com)
gmx rdf -s md.tpr -f "${TRAJ_FILE}" -n index.ndx -o rdf_pep_ache.xvg \
       -ref AChE -sel Peptide -selrpos mol_com -seltype mol_com -bin 0.02

# 自动读取轨迹总时间(ps)，自适应计算四等份区段，无视测试轨迹或正式100ns/1000ns差异
END_PS=$(gmx check -f "${TRAJ_FILE}" 2>&1 | grep -iE "Last frame|Step" | tail -n 1 | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+(\.[0-9]+)?$/) last=$i; print last}' || true)
if [ -z "${END_PS:-}" ] || [ "${END_PS}" = "0" ]; then
    END_PS=20
fi
STEP_PS=$(awk "BEGIN {print ${END_PS} / 4.0}")
echo ">> 自动读取轨迹长度: 0 - ${END_PS} ps, 四等份切分步长: ${STEP_PS} ps ..."

for q in 1 2 3 4; do
    t0=$(awk "BEGIN {print ($q - 1) * $STEP_PS}")
    t1=$(awk "BEGIN {print $q * $STEP_PS}")
    gmx rdf -s md.tpr -f "${TRAJ_FILE}" -n index.ndx \
        -o "rdf_pep_ache_q${q}.xvg" -ref AChE -sel Peptide \
        -selrpos mol_com -seltype mol_com -b "${t0}" -e "${t1}" -bin 0.02 || true
done

echo "RDF 输出: rdf_pep_ache.xvg 及 rdf_pep_ache_q1..4.xvg"
