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

$ErrorActionPreference = "Continue"

if ($OnlyMD) {
    $env:ONLY_MD = "1"
    Write-Host ">> [ONLY_MD Mode] Directly starting Production MD from Step 11/11..." -ForegroundColor Yellow
    & bash "./run_all.sh" $System "--only-md"
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} elseif ($Testing -or $env:TESTING -eq "1") {
    $env:ONLY_MD = "0"
    $env:TESTING = "1"
    Write-Host ">> [5000-Step Validation Mode] Running FULL Workflow from Step 1 (pdb2gmx) to Step 11 (MD 5000 steps) (mdp/test)" -ForegroundColor Yellow
    & bash "./run_all.sh" $System "--testing"
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
} else {
    $env:ONLY_MD = "0"
    $env:TESTING = "0"
    Write-Host ">> [100 ns Formal Simulation Mode] Running FULL 100 ns Production MD Workflow (mdp/100ns)" -ForegroundColor Cyan
    & bash "./run_all.sh" $System
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}
