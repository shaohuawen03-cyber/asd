#!/usr/bin/env bash
# ============================================================
# AChE (4ey6) - Aβ肽 对接复合物 分子动力学与论文分析 一键主流水线 (Bash / WSL 版)
#
# 用法:
#   ./run_pipeline_all.sh alllhrc              # 完整运行 100ns 模拟与全项分析
#   ONLY_MD=1 ./run_pipeline_all.sh alllhrc    # 跳过前10步平衡，直接跑 100ns MD 并执行分析
#   ONLY_ANALYSIS=1 ./run_pipeline_all.sh alllhrc # 已有轨迹，直接执行分析与绘图
# ============================================================
set -eu

SYS="${1:?用法: ./run_pipeline_all.sh <系统前缀, 如 alllhrc>}"

echo "===================================================================="
echo " 乙酰胆碱酯酶(4ey6)-Aβ肽 对接复合物 100 ns 分子动力学与全套分析 主流水线"
echo " 体系: ${SYS}"
echo "===================================================================="

if [ "${ONLY_ANALYSIS:-0}" != "1" ]; then
    echo ""
    echo ">>> [阶段一] 启动 100 ns 正式分子动力学模拟 (run_all.sh) ..."
    bash "./run_all.sh" "${SYS}"
else
    echo ""
    echo ">>> [阶段一] 仅分析模式：跳过动力学模拟阶段，直接从现有轨迹执行分析 ..."
fi

echo ""
echo ">>> [阶段二 & 阶段三] 运行对标论文多层次分析与批量生成 SVG/PNG/PDF 出版图表 (run_analysis.sh) ..."
bash "./run_analysis.sh" "${SYS}"

echo ""
echo "===================================================================="
echo " 大功告成！体系 ${SYS} 的 100 ns 动力学模拟及全套论文图表已全部顺利完成！"
echo " >> 出版级矢量/位图目录: md_${SYS}/figures/"
echo "    - fig0_summary_all.{svg,png,pdf} (A-F 2x3 综合主图)"
echo "    - fig1_rmsd_rmsf 到 fig6_bridging_waters 全套单图"
echo "    - summary_metrics.csv & wide.csv 统计指标表"
echo "===================================================================="
