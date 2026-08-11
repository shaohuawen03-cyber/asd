# ============================================================
# AChE - A-beta Complex MD & Analysis End-to-End Pipeline (PowerShell)
#
# Usage:  .\run_pipeline_all.ps1 -System alllhrc
#         .\run_pipeline_all.ps1 -System alllhrc -Testing
#         .\run_pipeline_all.ps1 -System alllhrc -OnlyMD
#         .\run_pipeline_all.ps1 -System alllhrc -OnlyAnalysis
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System,

    [switch]$Testing,
    [switch]$OnlyMD,
    [switch]$OnlyAnalysis
)

[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::InputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " 乙酰胆碱酯酶(4ey6) - Aβ肽 对接复合物 分子动力学与全套分析 主流水线" -ForegroundColor Green
Write-Host " 体系: $System" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

if (-not $OnlyAnalysis) {
    Write-Host ""
    Write-Host "=== 阶段一: 启动分子动力学模拟阶段 (run_all.ps1) ===" -ForegroundColor Cyan
    if ($Testing) {
        if ($OnlyMD) {
            & .\run_all.ps1 -System $System -Testing -OnlyMD 2>&1 | ForEach-Object { "$_" }
        } else {
            & .\run_all.ps1 -System $System -Testing 2>&1 | ForEach-Object { "$_" }
        }
    } else {
        if ($OnlyMD) {
            & .\run_all.ps1 -System $System -OnlyMD 2>&1 | ForEach-Object { "$_" }
        } else {
            & .\run_all.ps1 -System $System 2>&1 | ForEach-Object { "$_" }
        }
    }
    if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
        Write-Host "ERROR: run_all.ps1 failed with exit code $LASTEXITCODE" -ForegroundColor Red
        exit $LASTEXITCODE
    }
} else {
    Write-Host ""
    Write-Host "=== 阶段一: 仅分析模式 - 跳过动力学模拟阶段 ===" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "=== 阶段二 & 阶段三: 执行论文对标多层次分析与绘制出版级图表 (run_analysis.ps1) ===" -ForegroundColor Cyan
if ($Testing) {
    & .\run_analysis.ps1 -System $System -Testing 2>&1 | ForEach-Object { "$_" }
} else {
    & .\run_analysis.ps1 -System $System 2>&1 | ForEach-Object { "$_" }
}
if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
    Write-Host "ERROR: run_analysis.ps1 failed with exit code $LASTEXITCODE" -ForegroundColor Red
    exit $LASTEXITCODE
}

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " 大功告成！体系 $System MD 模拟与全套图表生成均顺利通关！" -ForegroundColor Green
Write-Host " 出版级矢量/位图保存目录: md_$System/figures/" -ForegroundColor Green
Write-Host "    - fig0_summary_all.{svg,png,pdf} (A-F 2x3 经典六格组合大总图)" -ForegroundColor Green
Write-Host "    - fig1_rmsd_rmsf 到 fig6_bridging_waters 全套独立指标单图" -ForegroundColor Green
Write-Host "    - summary_metrics.csv 及 wide.csv 统计指标表" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
