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

$ErrorActionPreference = "Stop"

if ($Testing -or $env:TESTING -eq "1") {
    Write-Host ">> [5000-Step Validation Mode] Running 5000-Step Testing Workflow (mdp/test)" -ForegroundColor Yellow
    $env:TESTING = "1"
} else {
    Write-Host ">> [100 ns Formal Simulation Mode] Running 100 ns Production MD Workflow (mdp/100ns)" -ForegroundColor Cyan
    $env:TESTING = "0"
}

if ($OnlyMD) {
    $env:ONLY_MD = "1"
    Write-Host ">> [ONLY_MD Mode] Directly starting Production MD from Step 11/11..." -ForegroundColor Yellow
    & bash "./run_all.sh" $System "--only-md"
} else {
    $env:ONLY_MD = "0"
    & bash "./run_all.sh" $System
}
