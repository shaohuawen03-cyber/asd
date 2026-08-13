# UTF-8. Console messages are English to avoid GBK mojibake on Windows.
# Rebuild COMPLEX DSSP + peptide RMSD diagnosis + replot only.
param(
    [Parameter(Mandatory=$true)]
    [string]$System
)

try { chcp 65001 | Out-Null } catch {}
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

$PY = "python"
if (Get-Command "python.exe" -ErrorAction SilentlyContinue) { $PY = "python.exe" }

$WorkDir = "..\md_$System"
if (-not (Test-Path $WorkDir)) {
    Write-Host "ERROR: missing $WorkDir" -ForegroundColor Red
    exit 1
}
if (-not (Test-Path "$WorkDir\md_0_1.xtc") -and -not (Test-Path "$WorkDir\md.xtc") -and -not (Test-Path "$WorkDir\md_fit.xtc")) {
    Write-Host "ERROR: no trajectory in $WorkDir" -ForegroundColor Red
    exit 1
}

Write-Host "========== DSSP(complex) + RMSD diagnosis: $System ==========" -ForegroundColor Green
Write-Host "Does NOT rerun RDF/SASA/Hbond/Contacts. Does NOT stop mdrun." -ForegroundColor Yellow

Push-Location $WorkDir
try {
    if (-not (Test-Path "index.ndx")) {
        Write-Host "[0] building index.ndx ..."
        & bash "../scripts/analysis/0_make_index.sh" 2>&1 | ForEach-Object { "$_" }
    } else {
        Write-Host "[0] index.ndx exists, skip trjconv"
    }

    Write-Host "[1/4] COMPLEX DSSP then peptide DSSP ..."
    & bash "../scripts/analysis/4_secondary_structure.sh" 2>&1 | ForEach-Object { "$_" }

    Write-Host "[2/4] Peptide RMSD jump diagnosis (PBC vs real) ..."
    & $PY "..\scripts\analysis\diagnose_peptide_rmsd.py" -d . 2>&1 | ForEach-Object { "$_" }

    Write-Host "[3/4] Optional phase table ..."
    & $PY "..\scripts\analysis\analyze_peptide_phases.py" -d . 2>&1 | ForEach-Object { "$_" }

    Write-Host "[4/4] Replot figures ..."
    & $PY "..\scripts\analysis\plot_all.py" --dir . --out ./figures 2>&1 | ForEach-Object { "$_" }
} finally {
    Pop-Location
}

Write-Host "========== done: md_$System\figures\ ==========" -ForegroundColor Green
Write-Host "Check:"
Write-Host "  ss_complex_summary.txt         (AChE+peptide DSSP, ~530 residues)"
Write-Host "  ss_pep_summary.txt             (peptide supplement only)"
Write-Host "  peptide_rmsd_jump_diagnosis.txt"
Write-Host "  fig4_secondary_structure.*     (complex helix/sheet, not 6-res coil)"
