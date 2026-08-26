# run_unified_replot.ps1
# One-click: unified replot of the 4 systems + the 3 compare_ache_vs_* folders.
# NO mdrun, NO MD re-run, NO md_* deletion. Safe while mdrun jobs are running.
#
# Usage (from F:\0wsh\asd\gromacs_md\scripts):
#   .\run_unified_replot.ps1

$ErrorActionPreference = "Stop"
try { chcp 65001 | Out-Null } catch {}
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8

$Here = $PSScriptRoot
Push-Location $Here
try {
    $py = Get-Command python -ErrorAction SilentlyContinue
    if (-not $py) { $py = Get-Command py -ErrorAction SilentlyContinue }
    if (-not $py) {
        Write-Host "ERROR: python not found on PATH (need python + matplotlib/numpy/pandas)" -ForegroundColor Red
        exit 1
    }
    Write-Host "python: $($py.Source)" -ForegroundColor Cyan
    & $py.Source "$Here\unified_replot_and_compare.py"
    exit $LASTEXITCODE
} finally {
    Pop-Location
}
