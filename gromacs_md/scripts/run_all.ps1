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
    Write-Host ">> [ONLY_MD Mode] Directly starting Production MD from Step 11/11..." -ForegroundColor Yellow
    & bash "./run_all.sh" $System "--only-md" 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} elseif ($Testing) {
    $env:ONLY_MD = "0"
    $env:TESTING = "1"
    Write-Host ">> [5000-Step Validation Mode] Running FULL Workflow from Step 1 (pdb2gmx) to Step 11 (MD 5000 steps) (mdp/test)" -ForegroundColor Yellow
    & bash "./run_all.sh" $System "--testing" 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} else {
    $env:ONLY_MD = "0"
    $env:TESTING = "0"
    $env:SIM_MODE = "100ns"
    Remove-Item env:TESTING -ErrorAction SilentlyContinue
    Remove-Item env:ONLY_MD -ErrorAction SilentlyContinue
    Write-Host ">> [100 ns Formal Simulation Mode] Running FULL 100 ns Production MD Workflow (mdp/100ns)" -ForegroundColor Cyan
    & bash "./run_all.sh" $System 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}
