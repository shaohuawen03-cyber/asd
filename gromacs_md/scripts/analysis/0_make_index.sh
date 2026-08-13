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

if [ -z "${ACHE_SEL:-}" ] || [ -z "${PEP_SEL:-}" ]; then
    if [ "${N_PROT_RES}" = "537" ]; then
        # 你的 7 肽对接构象体系: AChE(1-530), Peptide(531-537)
        ACHE_SEL="ri 1-530"
        PEP_SEL="ri 531-537"
    elif [ "${N_PROT_RES}" = "579" ]; then
        # 论文 Aβ(1-42) 体系: AChE(1-537), Peptide(538-579)
        ACHE_SEL="ri 1-537"
        PEP_SEL="ri 538-579"
    else
        # 自适应后退: 默认最后 7 个残基为小肽
        ACHE_SEL="ri 1-$((N_PROT_RES - 7))"
        PEP_SEL="ri $((N_PROT_RES - 6))-${N_PROT_RES}"
    fi
fi

# 自动检测该体系 tpr 中原有默认组的最大编号(如 16 或 18 等)
LAST_IDX=$(echo "q" | gmx make_ndx -f "${TPR_FILE}" -o /dev/null 2>&1 | grep -E "^ *Group +[0-9]+" | tail -n 1 | awk '{print $2}' || true)
LAST_IDX="${LAST_IDX:-16}"

G1=$((LAST_IDX + 1))
G2=$((LAST_IDX + 2))
G3=$((LAST_IDX + 3))
G4=$((LAST_IDX + 4))
G5=$((LAST_IDX + 5))

# 检查体系中是否含有第二条肽链或小肽残基 (识别是否为单独 AChE 单体对照组)
HAS_PEP=1
PEP_CNT=$(echo "${PEP_SEL}" | gmx make_ndx -f "${TPR_FILE}" -o /dev/null 2>&1 | grep -E "Found [1-9][0-9]* atoms" | wc -l || echo "0")
if [ "${PEP_CNT}" -eq "0" ]; then
    HAS_PEP=0
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
    echo ">> [单体蛋白模式] 未检测到第二条肽链，自动配置为单独 AChE 单体 (Protein alone) 对照组索引..."
    gmx make_ndx -f "${TPR_FILE}" -o index.ndx << EOF
${ACHE_SEL}
name ${G1} AChE
${G1} & 4
name ${G2} AChE_Backbone
q
EOF
    echo ">> 已成功生成单体 index.ndx, 包含组: AChE, AChE_Backbone。"
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
