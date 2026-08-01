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

# 自动将 Windows Anaconda / Miniconda 目录 (含 WSL /mnt/f 映射) 放入 PATH
for ana_dir in "/mnt/f/anaconda" "/mnt/f/Anaconda" "/f/anaconda" "/f/Anaconda" "/mnt/c/anaconda3" "/c/anaconda3" "$HOME/anaconda3" "$HOME/Anaconda3"; do
    if [ -d "$ana_dir" ]; then
        export PATH="$ana_dir:$ana_dir/Scripts:$ana_dir/Library/bin:$PATH"
        break
    fi
done

# 自动扫描并锁定已装有 MDAnalysis 的 Python 解释器 (兼容 WSL 下的 Windows /mnt/f/anaconda/python.exe)
PY_FOUND=""
if [ -n "${PYTHON:-}" ]; then
    PY_FOUND="${PYTHON}"
else
    for cand in \
        "/mnt/f/anaconda/python.exe" \
        "/mnt/f/Anaconda/python.exe" \
        "/mnt/c/anaconda3/python.exe" \
        "/f/anaconda/python.exe" \
        "/f/Anaconda/python.exe" \
        "/c/anaconda3/python.exe" \
        "python.exe" \
        "python" \
        "python3"; do
        if command -v "$cand" >/dev/null 2>&1 || [ -x "$cand" ]; then
            if "$cand" -c "import MDAnalysis" >/dev/null 2>&1; then
                PY_FOUND="$cand"
                break
            fi
        fi
    done
fi
PY="${PY_FOUND:-python}"
export PY
echo ">> [Python 环境] 锁定已装有 MDAnalysis 库的解释器: ${PY}"

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
echo "[0/8] 建立分析分组索引 (0_make_index.sh) ..."
bash "${SCRIPTS_DIR}/0_make_index.sh"

# 1. 骨架 RMSD / RMSF (论文 3.1 图1)
echo "[1/8] 计算骨架 RMSD / RMSF (1_rmsd_rmsf.sh) ..."
bash "${SCRIPTS_DIR}/1_rmsd_rmsf.sh"

# 2. 径向分布函数 RDF (论文 3.1 图2)
echo "[2/8] 计算 Aβ 围绕 AChE 的径向分布函数 RDF (2_rdf.sh) ..."
bash "${SCRIPTS_DIR}/2_rdf.sh"

# 3. 溶剂可及表面积 SASA (论文 3.2 图3)
echo "[3/8] 计算复合物 SASA (3_sasa.sh) ..."
bash "${SCRIPTS_DIR}/3_sasa.sh"

# 4. 二级结构分析 DSSP (论文 3.2 图4)
echo "[4/8] 计算 Aβ 肽二级结构演变 (4_secondary_structure.sh) ..."
bash "${SCRIPTS_DIR}/4_secondary_structure.sh"

# 5. 氢键数量统计 (论文 3.3)
echo "[5/8] 统计间/内氢键分布 (5_hbond.sh) ..."
bash "${SCRIPTS_DIR}/5_hbond.sh"

to_py_path() {
    local p="$1"
    if command -v wslpath >/dev/null 2>&1 && [[ "${PY}" == *".exe"* || "${PY}" == *"/mnt/"* || "${PY}" == *":"* ]]; then
        wslpath -w "$p" 2>/dev/null || echo "$p"
    else
        echo "$p"
    fi
}

# 6. 非天然接触统计 (论文 3.3 图5/表1)
echo "[6/8] 统计天然与非天然相互作用接触 (contacts.py) ..."
"${PY}" "$(to_py_path "${SCRIPTS_DIR}/contacts.py")" -t md.tpr -f md.xtc

# 7. 水介导桥连相互作用 (论文 3.4 图6/表2)
echo "[7/8] 统计水介导桥连相互作用 (bridging_waters.py) ..."
"${PY}" "$(to_py_path "${SCRIPTS_DIR}/bridging_waters.py")" -t md.tpr -f md.xtc

# 8. 批量生成矢量/位图出版图与统计表 (plot_all.py)
echo "[8/8] 批量绘制论文出版级图表 (SVG / PNG / PDF) ..."
"${PY}" "$(to_py_path "${SCRIPTS_DIR}/plot_all.py")" --dir . --out ./figures

echo "========== 体系 ${SYS} 分析流程全部完成！ =========="
echo "生成分析图表一览 (保存在 ./figures/ 下):"
echo "  - fig1_rmsd_rmsf.{svg,png,pdf}            => 图 1 (RMSD / RMSF)"
echo "  - fig2_rdf.{svg,png,pdf}                  => 图 2 (径向分布函数 RDF)"
echo "  - fig3_sasa.{svg,png,pdf}                 => 图 3 (溶剂可及表面积 SASA)"
echo "  - fig4_secondary_structure.{svg,png,pdf}  => 图 4 (二级结构倾向)"
echo "  - fig5_contacts.{svg,png,pdf}             => 图 5 (非天然残基对接触)"
echo "  - fig6_bridging_waters.{svg,png,pdf}      => 图 6 (水介导桥连)"
echo "  - fig_hbonds.{svg,png,pdf}                => 论文 3.3 节 (氢键数量曲线)"
echo "  - fig0_summary_all.{svg,png,pdf}          => 综合 2x3 六格汇总对比主图"
echo "  - summary_metrics.csv & wide.csv          => 统计指标汇总表"
