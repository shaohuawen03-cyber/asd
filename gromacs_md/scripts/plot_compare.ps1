# ============================================================
# AChE Monomer vs Complex Comparison Plotting Suite (PowerShell)
#
# Usage:  .\plot_compare.ps1 -Protein ache -Complex alllhrc
# ============================================================
param(
    [string]$Protein = "ache",
    [string]$Complex = "alllhrc",
    [string]$OutDir = "..\compare_figures"
)

[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::InputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

$PY = "python"
if (Get-Command "python.exe" -ErrorAction SilentlyContinue) {
    $PY = "python.exe"
}

$ScriptsDir = (Resolve-Path ".\analysis").Path

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " AChE Monomer ($Protein) vs Complex ($Complex) Comparative Plotting" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

& $PY "$ScriptsDir\plot_compare_systems.py" --protein "..\md_$Protein" --complex "..\md_$Complex" --out $OutDir 2>&1 | ForEach-Object { "$_" }
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " SUCCESS! Comparative figures saved in $OutDir\" -ForegroundColor Green
Write-Host "   - protein_vs_complex.{svg,png,pdf} (2x3 Combined Master Plot)" -ForegroundColor Green
Write-Host "   - compare_summary_metrics.csv" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
