# ============================================================
# AChE - A-beta Complex 100 ns Formal Production MD & Analysis Suite (PowerShell)
#
# Usage:  .\run_100ns_formal.ps1 -System alllhrc
#         .\run_100ns_formal.ps1 -System fllhttr
#         .\run_100ns_formal.ps1 -System ylsllqr
#         .\run_100ns_formal.ps1 -System ache
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System
)

[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::InputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

# 强制清除之前任何测试会话残留的环境变量，1000% 锁定为 100 ns 正式产物生产模式
$env:TESTING = "0"
$env:ONLY_MD = "0"
$env:SIM_MODE = "100ns"
Remove-Item env:TESTING -ErrorAction SilentlyContinue
Remove-Item env:ONLY_MD -ErrorAction SilentlyContinue

$WorkDir = "..\md_$System"

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " AChE(4ey6) - A-beta Complex 100 ns Formal Production MD & Analysis Pipeline" -ForegroundColor Green
Write-Host " System: $System" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

Write-Host ""
Write-Host ">>> [Step 0/3] Deleting any previous test results in $WorkDir to ensure a 100% clean start..." -ForegroundColor Yellow
if (Test-Path $WorkDir) {
    Remove-Item -Path $WorkDir -Recurse -Force -ErrorAction SilentlyContinue
    Write-Host "  [DELETED PREVIOUS DIRECTORY] $WorkDir" -ForegroundColor Yellow
}

Write-Host ""
Write-Host ">>> [Step 1/3] Starting 100 ns Formal Production MD Simulation (mdp/100ns) ..." -ForegroundColor Cyan
& .\run_all.ps1 -System $System
if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
    Write-Host "ERROR: run_all.ps1 failed with exit code $LASTEXITCODE" -ForegroundColor Red
    exit $LASTEXITCODE
}

Write-Host ""
Write-Host ">>> [Step 2/3 & 3/3] Running Trajectory Analysis & Generating Publication SVG/PNG/PDF Figures ..." -ForegroundColor Cyan
& .\run_analysis.ps1 -System $System
if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
    Write-Host "ERROR: run_analysis.ps1 failed with exit code $LASTEXITCODE" -ForegroundColor Red
    exit $LASTEXITCODE
}

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " SUCCESS! System $System 100 ns Production MD simulation and analysis completed!" -ForegroundColor Green
Write-Host " Publication Figures Directory: md_$System/figures/" -ForegroundColor Green
Write-Host "    - fig0_summary_all.{svg,png,pdf} (A-F 2x3 Combined Master Plot)" -ForegroundColor Green
Write-Host "    - fig1_rmsd_rmsf to fig6_bridging_waters individual plots" -ForegroundColor Green
Write-Host "    - summary_metrics.csv and wide.csv statistical summary tables" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
