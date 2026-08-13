# Shared replot / re-analysis for ONE finished v1.0 system.
# Does NOT start mdrun. Does NOT delete md_*.
#
# Usage (normally called by replot_<system>.ps1):
#   .\replot_common.ps1 -System alllhrc
#   .\replot_common.ps1 -System alllhrc -OnlyPlot
#   .\replot_common.ps1 -System alllhrc -Full
#
# Default (no switch):
#   * if RMSD + complex DSSP already exist -> OnlyPlot (fast; use after a plot crash)
#   * if gyrate_complex.xvg is missing     -> run 6_rg.sh then plot
#   * otherwise                            -> full analysis (PBC / DSSP / Rg / H-bonds / ...)
param(
    [Parameter(Mandatory = $true)]
    [string]$System,

    [switch]$OnlyPlot,
    [switch]$Full
)

try { chcp 65001 | Out-Null } catch {}
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

$Here = $PSScriptRoot
$Work = Join-Path (Split-Path $Here -Parent) "md_$System"

Write-Host "============================================================" -ForegroundColor Cyan
Write-Host " REPLOT $System   (v1.0 trajectory + current analysis code)" -ForegroundColor Cyan
Write-Host " Does NOT start mdrun. Does NOT delete md_*." -ForegroundColor Yellow
Write-Host "============================================================" -ForegroundColor Cyan

if (-not (Test-Path $Work)) {
    Write-Host "ERROR: missing $Work" -ForegroundColor Red
    exit 1
}

$HasXtc = (Test-Path (Join-Path $Work "md.xtc")) -or (Test-Path (Join-Path $Work "md_0_1.xtc"))
if (-not $HasXtc) {
    Write-Host "ERROR: no md.xtc / md_0_1.xtc in $Work (MD not finished or not started)" -ForegroundColor Red
    exit 1
}

$HasGro = (Test-Path (Join-Path $Work "md.gro")) -or (Test-Path (Join-Path $Work "md_0_1.gro"))
if (-not $HasGro) {
    Write-Host "SKIP: xtc exists but no md.gro -- mdrun is probably still running." -ForegroundColor Yellow
    Write-Host "      Wait until md.gro appears, then re-run this script." -ForegroundColor Yellow
    exit 2
}

$HasRmsd = Test-Path (Join-Path $Work "rmsd_complex_bb.xvg")
$HasDssp = Test-Path (Join-Path $Work "ss_complex_frac.xvg")
$HasRg   = Test-Path (Join-Path $Work "gyrate_complex.xvg")

$Mode = "full"
if ($OnlyPlot -and $Full) {
    Write-Host "ERROR: do not pass -OnlyPlot and -Full together" -ForegroundColor Red
    exit 1
} elseif ($OnlyPlot) {
    $Mode = "plot"
} elseif ($Full) {
    $Mode = "full"
} elseif ($HasRmsd -and $HasDssp) {
    $Mode = "plot"
    if (-not $HasRg) { $Mode = "rg+plot" }
    Write-Host "Found existing RMSD + complex DSSP -> $Mode (pass -Full to recompute all analysis)" -ForegroundColor Yellow
} else {
    $Mode = "full"
    Write-Host "Missing RMSD or complex DSSP -> full analysis" -ForegroundColor Yellow
}

$Ana = Join-Path $Here "run_analysis.ps1"
$code = 0

# run_analysis.ps1 uses ..\md_$System relative to CWD — always run from scripts/
Push-Location $Here
try {
    if ($Mode -eq "rg+plot") {
        Write-Host "[Rg] gyrate_complex.xvg missing; running 6_rg.sh only ..." -ForegroundColor Cyan
        Push-Location $Work
        try {
            & bash "../scripts/analysis/6_rg.sh" 2>&1 | ForEach-Object { "$_" }
        } finally {
            Pop-Location
        }
        $Mode = "plot"
    }

    if ($Mode -eq "plot") {
        Write-Host "[Plot] run_analysis.ps1 -System $System -OnlyPlot" -ForegroundColor Cyan
        & $Ana -System $System -OnlyPlot
    } else {
        Write-Host "[Full] run_analysis.ps1 -System $System" -ForegroundColor Cyan
        & $Ana -System $System
    }
    $code = $LASTEXITCODE
} finally {
    Pop-Location
}

if ($code -ne 0 -and $null -ne $code) {
    Write-Host "[FAIL] $System exit $code" -ForegroundColor Red
    exit $code
}

Write-Host "[OK] $System -> md_$System\figures\" -ForegroundColor Green
Write-Host "Check:"
Write-Host "  ss_complex_summary.txt     (~530 residues, 0-100 ns, helix ~30%+)"
Write-Host "  figures\fig4_secondary_structure.png   DSSP % lines + occupancy bars + heatmap + peptide"
Write-Host "  figures\fig_peptide_rmsd_rmsf.png      peptide RMSD + RMSF (not on overview)"
Write-Host "  figures\fig0_summary_all.png           A RMSD B AChE-RMSF C RDF D SASA E DSSP-bars F Rg G H-bonds H DSSP-lines"
exit 0
