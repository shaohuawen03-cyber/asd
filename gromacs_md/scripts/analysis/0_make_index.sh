#!/usr/bin/env bash
# ============================================================
# 建立用于分析的分组索引 (AChE / 肽 / 骨架 / 水)
#
# v2.7.3: 分组直接由 pdb2gmx 拓扑推导, 与 PDB 残基编号完全无关:
#   - topol_Protein_chain_A.itp 的 [ atoms ] 行数 = NA (AChE 链原子数)
#   - topol_Protein_chain_B.itp 的 [ atoms ] 行数 = NB (肽链原子数; 无肽时文件不存在)
#   - tpr 中蛋白原子按 itp 顺序连续: AChE = 原子 1..NA, 肽 = 原子 NA+1..NA+NB
# 这样 4-542+1-7 编号的 PDB 与 1-530+531-537 编号的 PDB 都能正确分组。
#
# 用法: 在产物工作目录中运行
# ============================================================
set -eu

if [ -n "${GMX:-}" ]; then
    gmx() { "${GMX}" "$@"; }
elif command -v gmx.exe >/dev/null 2>&1; then
    gmx() { gmx.exe "$@"; }
fi

TPR_FILE="md.tpr"
if [ ! -f "${TPR_FILE}" ]; then
    TPR_FILE="md_0_1.tpr"
    if [ ! -f "${TPR_FILE}" ]; then
        TPR_FILE="neutral.gro"
    fi
fi

# 自动检测体系的蛋白质总残基数 (由于 binary tpr 不保存字母 chainID, 这里用残基号划分最稳定)
N_PROT_RES=$(echo "q" | gmx make_ndx -f "${TPR_FILE}" -o /dev/null 2>&1 | grep "Protein residues" | grep -oE "[0-9]+" | head -n 1 || true)
N_PROT_RES="${N_PROT_RES:-537}"

# ---------- 从 pdb2gmx 拓扑推导链组成 (编号无关, v2.7.3) ----------
count_itp_atoms() {
    local f="$1"
    [ -f "$f" ] || { echo 0; return; }
    awk '/^[[:space:]]*\[[[:space:]]*atoms[[:space:]]*\]/{f=1; next}
         /^[[:space:]]*\[/{f=0}
         f && $1 ~ /^[0-9]+$/ {c++}
         END{print c+0}' "$f"
}

NA=$(count_itp_atoms "topol_Protein_chain_A.itp")
NB=$(count_itp_atoms "topol_Protein_chain_B.itp")
NP=$((NA + NB))

HAS_PEP=1
if [ ! -f "topol_Protein_chain_B.itp" ] || [ "${NB}" -lt 2 ]; then
    HAS_PEP=0
fi

if [ "${NP}" -lt 1 ]; then
    # 极少数没有 topol itp 的环境: 退回老逻辑 (残基编号, 仅最后手段)
    echo "!! 找不到 topol_Protein_chain_A.itp，回退到残基编号猜测 (可能不准确) ..."
    HAS_PEP=1
    case "${N_PROT_RES}" in
        530) HAS_PEP=0; ACHE_SEL="ri 1-530" ;;
        537) ACHE_SEL="ri 1-530"; PEP_SEL="ri 531-537" ;;
        579) ACHE_SEL="ri 1-537"; PEP_SEL="ri 538-579" ;;
        *)   ACHE_SEL="ri 1-$((N_PROT_RES - 7))"; PEP_SEL="ri $((N_PROT_RES - 6))-${N_PROT_RES}" ;;
    esac
else
    # 标准路径: 原子区间分组, 编号无关
    ACHE_SEL="a 1-${NA}"
    if [ "${HAS_PEP}" = "1" ]; then
        PEP_SEL="a $((NA + 1))-${NP}"
    fi
    echo ">> [topol] chain A = ${NA} atoms; chain B = ${NB} atoms; N_PROT_RES(提示) = ${N_PROT_RES}"
fi

# 自动检测该体系 tpr 中原有默认组的最大编号(如 16 或 18 等)
LAST_IDX=$(echo "q" | gmx make_ndx -f "${TPR_FILE}" -o /dev/null 2>&1 | grep -E "^ *Group +[0-9]+" | tail -n 1 | awk '{print $2}' || true)
LAST_IDX="${LAST_IDX:-16}"

G1=$((LAST_IDX + 1))
G2=$((LAST_IDX + 2))
G3=$((LAST_IDX + 3))
G4=$((LAST_IDX + 4))
G5=$((LAST_IDX + 5))

# 复合物分支再做一次校验: 若 PEP_SEL 实际选不到任何原子, 降级为单体模式 (防编号偏差)
if [ "${HAS_PEP}" = "1" ]; then
    PEP_CNT=$(echo "${PEP_SEL}" | gmx make_ndx -f "${TPR_FILE}" -o /dev/null 2>&1 | grep -E "Found [1-9][0-9]* atoms" | wc -l || echo "0")
    if [ "${PEP_CNT}" -eq "0" ]; then
        HAS_PEP=0
    fi
fi
if [ "${HAS_PEP}" = "1" ] && [ "${NB}" -ge 2 ]; then
    # 拓扑链 B 原子数应等于 Peptide 组原子数; 不一致说明原子区间推导有误
    PEP_ATOMS=$(echo "${PEP_SEL}" | gmx make_ndx -f "${TPR_FILE}" -o /dev/null 2>&1 | grep -oE "Found [0-9]+ atoms" | grep -oE "[0-9]+" | head -n 1 || echo "0")
    if [ "${PEP_ATOMS:-0}" -ne "${NB}" ]; then
        echo "!! [警告] Peptide 组原子数(${PEP_ATOMS})与拓扑 chain B(${NB})不一致，但仍按拓扑分组继续。"
    fi
fi

if [ "${HAS_PEP}" = "1" ]; then
    echo ">> [复合物模式] 分配新增分组编号: ${G1}(AChE), ${G2}(Peptide), ${G3}(AChE_Backbone), ${G4}(Peptide_Backbone), ${G5}(Peptide_Calpha) ..."
    gmx make_ndx -f "${TPR_FILE}" -o index.ndx << EOF
${ACHE_SEL}
name ${G1} AChE
${PEP_SEL}
name ${G2} Peptide
${G1} & 4
name ${G3} AChE_Backbone
${G2} & 4
name ${G4} Peptide_Backbone
${G2} & 3
name ${G5} Peptide_Calpha
q
EOF
    echo ">> 已成功生成复合物 index.ndx, 包含组: AChE, Peptide, AChE_Backbone, Peptide_Backbone, Peptide_Calpha。"
else
    echo ">> [单体蛋白模式] 无肽 (单独 AChE 单体对照, 拓扑中无 chain B)，自动配置单体 index.ndx ..."
    gmx make_ndx -f "${TPR_FILE}" -o index.ndx << EOF
${ACHE_SEL}
name ${G1} AChE
${G1} & 4
name ${G2} AChE_Backbone
q
EOF
    echo ">> 已成功生成单体 index.ndx, 包含组: AChE, AChE_Backbone (无 Peptide 组)。"
    echo ">> [单体蛋白模式] 清理早期误算残留的肽相关分析产物 (单体对照不参与 AChE-肽氢键/RDF/肽 RMSD 等对比) ..."
    rm -f hbond_ache_pep.xvg hbond_pep_intra.xvg \
          rdf_pep_ache.xvg rdf_pep_ache_q1.xvg rdf_pep_ache_q2.xvg \
          rdf_pep_ache_q3.xvg rdf_pep_ache_q4.xvg \
          rmsd_pep_bb.xvg rmsf_pep_bb.xvg rmsd_pep_on_ache.xvg \
          sasa_pep.xvg \
          ss_pep.dat ss_pep_num.xvg ss_pep_frac.xvg ss_pep_bins.dat \
          ss_pep_perres.dat ss_pep_rama.dat ss_pep_summary.txt \
          inter_contacts.csv intra_contacts.csv frequent_contacts.tsv \
          bridging_per_residue.csv \
          peptide_rmsd_phases.dat peptide_rmsd_phases.txt \
          peptide_rmsd_jump_diagnosis.txt 2>/dev/null || true
    echo ">> [单体蛋白模式] 肽相关产物清理完成。"
fi

# ----- 去除周期性边界条件(PBC)、消除多链跨界拆分并叠合主干旋转平移 -----
RAW_XTC=""
if [ -f "md.xtc" ]; then
    RAW_XTC="md.xtc"
elif [ -f "md_0_1.xtc" ]; then
    RAW_XTC="md_0_1.xtc"
fi

if [ -n "${RAW_XTC}" ]; then
    if [ ! -f "md_fit.xtc" ] || [ "${RAW_XTC}" -nt "md_fit.xtc" ]; then
        echo ">> [1/4 恢复多肽与蛋白分子完整性] 正在执行 gmx trjconv -pbc whole ..."
        gmx trjconv -s "${TPR_FILE}" -f "${RAW_XTC}" -n index.ndx -o md_whole.xtc -pbc whole << EOF
0
EOF
        echo ">> [2/4 消除跨界折回突跳] 正在执行 gmx trjconv -pbc nojump 保证两链连续不跳变 ..."
        gmx trjconv -s "${TPR_FILE}" -f md_whole.xtc -n index.ndx -o md_nojump.xtc -pbc nojump << EOF
0
EOF
        echo ">> [3/4 紧凑居中与去PBC] 正在执行 gmx trjconv -center -pbc mol -ur compact 生成 md_center.xtc ..."
        gmx trjconv -s "${TPR_FILE}" -f md_nojump.xtc -n index.ndx -o md_center.xtc -center -pbc mol -ur compact << EOF
1
0
EOF
        echo ">> [4/4 旋转平移叠合] 正在执行 gmx trjconv -fit rot+trans 消除主干整体漂移生成 md_fit.xtc ..."
        gmx trjconv -s "${TPR_FILE}" -f md_center.xtc -n index.ndx -o md_fit.xtc -fit rot+trans << EOF
4
0
EOF
        rm -f md_whole.xtc md_nojump.xtc md_center.xtc 2>/dev/null || true
        echo ">> [OK] 已成功生成彻底去 PBC、无跨越突跳、紧凑居中且主干对齐的纯净轨迹: md_fit.xtc ！"
    fi
fi
