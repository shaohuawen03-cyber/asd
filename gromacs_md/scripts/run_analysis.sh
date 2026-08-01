#!/usr/bin/env bash
# ============================================================
# AChE (4ey6) - beta-淀粉样肽(Aβ)复合物 分子动力学模拟 分析全流程
# 复现论文《乙酰胆碱酯酶-β-淀粉样肽复合物的分子动力学模拟》第 3 部分
#
# 用法:  ./run_analysis.sh <系统前缀>
#        例:  TESTING=1 ./run_analysis.sh alllhrc
#              ./run_analysis.sh alllhrc
# ============================================================
set -eu

# 自动检测 gmx 或 gmx.exe
if [ -n "${GMX:-}" ]; then
    GMX_CMD="${GMX}"
elif command -v gmx.exe >/dev/null 2>&1; then
    GMX_CMD="gmx.exe"
elif command -v gmx >/dev/null 2>&1; then
    GMX_CMD="gmx"
else
    GMX_CMD="gmx"
fi
export GMX="${GMX_CMD}"

SYS="${1:?用法: ./run_analysis.sh <前缀, 如 alllhrc>}"
WORK="../md_${SYS}"

if [ ! -d "${WORK}" ] || [ ! -f "${WORK}/md.xtc" ]; then
    echo "!!! 错误: 在 ${WORK} 目录下未找到产物轨迹文件 md.xtc！"
    echo "!!! 请先执行 ./run_all.sh ${SYS} 完成 MD 模拟。"
    exit 1
fi

export TESTING="${TESTING:-0}"
if [ "${TESTING}" = "1" ]; then
    echo ">> [分析测试模式] 针对短时间测试轨迹自动自适应参数"
else
    echo ">> [分析正式模式] 针对 1000 ns 生产轨迹执行完整计算"
fi

cd "${WORK}"
# 取得相对脚本绝对路径
SCRIPTS_DIR="$(cd ../scripts/analysis && pwd)"

echo "========== 开始对体系 ${SYS} 执行论文结果分析 ========"

# 0. 建立分析索引
echo "[0/7] 建立分析分组索引 (0_make_index.sh) ..."
bash "${SCRIPTS_DIR}/0_make_index.sh"

# 1. 骨架 RMSD / RMSF (论文 3.1 图1)
echo "[1/7] 计算骨架 RMSD / RMSF (1_rmsd_rmsf.sh) ..."
bash "${SCRIPTS_DIR}/1_rmsd_rmsf.sh"

# 2. 径向分布函数 RDF (论文 3.1 图2)
echo "[2/7] 计算 Aβ 围绕 AChE 的径向分布函数 RDF (2_rdf.sh) ..."
bash "${SCRIPTS_DIR}/2_rdf.sh"

# 3. 溶剂可及表面积 SASA (论文 3.2 图3)
echo "[3/7] 计算复合物 SASA (3_sasa.sh) ..."
bash "${SCRIPTS_DIR}/3_sasa.sh"

# 4. 二级结构分析 DSSP (论文 3.2 图4)
echo "[4/7] 计算 Aβ 肽二级结构演变 (4_secondary_structure.sh) ..."
bash "${SCRIPTS_DIR}/4_secondary_structure.sh"

# 5. 氢键数量统计 (论文 3.3)
echo "[5/7] 统计间/内氢键分布 (5_hbond.sh) ..."
bash "${SCRIPTS_DIR}/5_hbond.sh"

# 6. 非天然接触统计 (论文 3.3 图5/表1)
echo "[6/7] 统计天然与非天然相互作用接触 (contacts.py) ..."
if python3 -c "import MDAnalysis" >/dev/null 2>&1; then
    python3 "${SCRIPTS_DIR}/contacts.py" -t md.tpr -f md.xtc
else
    echo ">> [提示] 未检测到 Python MDAnalysis 库，跳过 contacts.py 计算。"
    echo ">>        如需执行，请通过命令: pip install -r ../scripts/analysis/requirements.txt 安装。"
fi

# 7. 水介导桥连相互作用 (论文 3.4 图6/表2)
echo "[7/7] 统计水介导桥连相互作用 (bridging_waters.py) ..."
if python3 -c "import MDAnalysis" >/dev/null 2>&1; then
    python3 "${SCRIPTS_DIR}/bridging_waters.py" -t md.tpr -f md.xtc
else
    echo ">> [提示] 未检测到 Python MDAnalysis 库，跳过 bridging_waters.py 计算。"
fi

echo "========== 体系 ${SYS} 分析流程全部完成！ =========="
echo "生成分析产物对照论文表一览:"
echo "  - rmsd_*_bb.xvg, rmsf_*_bb.xvg            => 图 1 (RMSD / RMSF)"
echo "  - rdf_pep_ache*.xvg                       => 图 2 (径向分布函数 RDF)"
echo "  - sasa_*.xvg                              => 图 3 (溶剂可及表面积 SASA)"
echo "  - ss_pep.xpm, ss_pep.sc, ss_pep_bins.dat  => 图 4 (二级结构倾向)"
echo "  - hbond_*.xvg                             => 论文 3.3 节 (氢键分布)"
echo "  - inter/intra_contacts.csv                => 图 5 / 表 1 (非天然接触)"
echo "  - bridging_per_residue.csv                => 图 6 / 表 2 (桥连水分析)"
