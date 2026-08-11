# ============================================================
# AChE (4ey6) - beta-amyloid peptide MD Workflow (PowerShell)
#
# Usage:  .\run_all.ps1 -System alllhrc -Testing
#         .\run_all.ps1 -System alllhrc -OnlyMD
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System,

    [switch]$Testing,
    [switch]$OnlyMD
)

[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::InputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

if ($OnlyMD) {
    $env:ONLY_MD = "1"
    Write-Host ">> [ONLY_MD Mode] 仅启动第 [11/11] 步产物动力学模拟 ..." -ForegroundColor Yellow
    & bash "./run_all.sh" $System "--only-md" 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} elseif ($Testing -or $env:TESTING -eq "1") {
    $env:ONLY_MD = "0"
    $env:TESTING = "1"
    Write-Host ">> [5,000步验证模式] 完整运行所有阶段 (准备结构->pdb2gmx->三斜盒子->加水->生理盐->索引->EM->退火->NPT->5000步MD)" -ForegroundColor Yellow
    & bash "./run_all.sh" $System "--testing" 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} else {
    $env:ONLY_MD = "0"
    $env:TESTING = "0"
    Write-Host ">> [100 ns 正式生产模拟模式] 完整运行 100 ns 生产级流程度 (mdp/100ns)" -ForegroundColor Cyan
    & bash "./run_all.sh" $System 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}
