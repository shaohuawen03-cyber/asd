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

$ErrorActionPreference = "Continue"

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " AChE(4ey6) - A-beta Complex MD & Analysis Master Pipeline" -ForegroundColor Green
Write-Host " System: $System" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

if (-not $OnlyAnalysis) {
    Write-Host ""
    Write-Host "=== Phase 1: Starting Molecular Dynamics Simulation (run_all.ps1) ===" -ForegroundColor Cyan
    if ($Testing) {
        if ($OnlyMD) {
            .\run_all.ps1 -System $System -Testing -OnlyMD
        } else {
            .\run_all.ps1 -System $System -Testing
        }
    } else {
        if ($OnlyMD) {
            .\run_all.ps1 -System $System -OnlyMD
        } else {
            .\run_all.ps1 -System $System
        }
    }
    if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
        Write-Host "ERROR: run_all.ps1 failed with exit code $LASTEXITCODE" -ForegroundColor Red
        exit $LASTEXITCODE
    }
} else {
    Write-Host ""
    Write-Host "=== Phase 1: OnlyAnalysis mode - Skipping MD simulation ===" -ForegroundColor Yellow
}

Write-Host ""
Write-Host "=== Phase 2 & 3: Running Trajectory Analysis & Generating Figures ===" -ForegroundColor Cyan
if ($Testing) {
    .\run_analysis.ps1 -System $System -Testing
} else {
    .\run_analysis.ps1 -System $System
}
if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
    Write-Host "ERROR: run_analysis.ps1 failed with exit code $LASTEXITCODE" -ForegroundColor Red
    exit $LASTEXITCODE
}

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " SUCCESS! System $System MD simulation and analysis plots completed!" -ForegroundColor Green
Write-Host " Publication Figures Directory: md_$System/figures/" -ForegroundColor Green
Write-Host "    - fig0_summary_all.{svg,png,pdf} (A-F 2x3 Combined Master Plot)" -ForegroundColor Green
Write-Host "    - fig1_rmsd_rmsf to fig6_bridging_waters individual plots" -ForegroundColor Green
Write-Host "    - summary_metrics.csv and wide.csv statistical summary tables" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
