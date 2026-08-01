# ============================================================
# AChE (4ey6) - beta-amyloid peptide MD Workflow (PowerShell)
#
# Usage:  .\run_all.ps1 -System alllhrc -Testing
#         .\run_all.ps1 -System alllhrc
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System,

    [switch]$Testing
)

$ErrorActionPreference = "Stop"

if ($Testing -or $env:TESTING -eq "1") {
    Write-Host ">> [TESTING MODE] 100-step short MD test workflow" -ForegroundColor Yellow
    $env:TESTING = "1"
} else {
    Write-Host ">> [PRODUCTION MODE] 1000 ns production MD workflow" -ForegroundColor Cyan
    $env:TESTING = "0"
}

& bash ".\run_all.sh" $System
