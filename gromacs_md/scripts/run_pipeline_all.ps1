# ============================================================
# AChE (4ey6) - Aβ肽 对接复合物 分子动力学与论文分析 一键主流水线 (PowerShell 原生版)
#
# 从头到尾一键全自动跑完:
#   1. 100 ns 正式动力学模拟 (三斜盒子优化、AMBER amber99sb-ildn/ff14SB 力场、0->300K 连续退火)
#   2. 7 大论文对标轨迹统计与计算 (RMSD/RMSF/RDF/SASA/DSSP/Hbonds/Contacts/Bridging Waters)
#   3. 批量生成发表级 SVG / PNG / PDF 图像与表 1 - 表 2 统计数据表
#
# 用法:
#   .\run_pipeline_all.ps1 -System alllhrc              # 从头开始完整运行 100ns 模拟与全项分析
#   .\run_pipeline_all.ps1 -System alllhrc -OnlyMD      # 跳过前10步平衡，直接跑 100ns MD 并执行分析
#   .\run_pipeline_all.ps1 -System alllhrc -OnlyAnalysis # 已有 md.xtc 轨迹，直接执行全流程分析绘图
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System,

    [switch]$OnlyMD,
    [switch]$OnlyAnalysis
)

$ErrorActionPreference = "Stop"

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " 乙酰胆碱酯酶(4ey6)-Aβ肽 对接复合物 100 ns 分子动力学与全套分析 主流水线" -ForegroundColor Green
Write-Host " 体系: $System" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

if (-not $OnlyAnalysis) {
    Write-Host "`n>>> [阶段一] 启动 100 ns 正式分子动力学模拟 (run_all.ps1) ..." -ForegroundColor Cyan
    if ($OnlyMD) {
        .\run_all.ps1 -System $System -OnlyMD
    } else {
        .\run_all.ps1 -System $System
    }
} else {
    Write-Host "`n>>> [阶段一] 仅分析模式：跳过动力学模拟阶段，直接从现有轨迹执行分析 ..." -ForegroundColor Yellow
}

Write-Host "`n>>> [阶段二 & 阶段三] 运行对标论文多层次分析与批量生成 SVG/PNG/PDF 出版图表 (run_analysis.ps1) ..." -ForegroundColor Cyan
.\run_analysis.ps1 -System $System

Write-Host "`n====================================================================" -ForegroundColor Green
Write-Host " 大功告成！体系 $System 的 100 ns 动力学模拟及全套论文图表已全部顺利完成！" -ForegroundColor Green
Write-Host " >> 出版级矢量/位图目录: md_$System\figures\" -ForegroundColor Green
Write-Host "    - fig0_summary_all.{svg,png,pdf} (A-F 2x3 综合主图)" -ForegroundColor Green
Write-Host "    - fig1_rmsd_rmsf 到 fig6_bridging_waters 全套单图" -ForegroundColor Green
Write-Host "    - summary_metrics.csv & wide.csv 统计指标表" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
