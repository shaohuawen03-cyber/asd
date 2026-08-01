# ============================================================
# AChE (4ey6) - beta-amyloid peptide MD Analysis Workflow (PowerShell)
# Reproducing Paper Section 3 Analyses
#
# Usage:  .\run_analysis.ps1 -System alllhrc -Testing
#         .\run_analysis.ps1 -System alllhrc
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System,

    [switch]$Testing
)

$ErrorActionPreference = "Stop"

$GMX = "gmx"
if (Get-Command "gmx.exe" -ErrorAction SilentlyContinue) {
    $GMX = "gmx.exe"
}

$PY = "python"
if (Get-Command "python.exe" -ErrorAction SilentlyContinue) {
    $PY = "python.exe"
}

$WorkDir = "..\md_$System"
if (-not (Test-Path "$WorkDir\md.xtc")) {
    Write-Host "ERROR: Trajectory file md.xtc not found in $WorkDir!" -ForegroundColor Red
    Write-Host "Please run MD simulation first." -ForegroundColor Red
    exit 1
}

if ($Testing -or $env:TESTING -eq "1") {
    Write-Host ">> [TESTING MODE] 100-step short trajectory parameters" -ForegroundColor Yellow
    $env:TESTING = "1"
} else {
    Write-Host ">> [PRODUCTION MODE] 1000 ns production trajectory analysis" -ForegroundColor Cyan
    $env:TESTING = "0"
}

Push-Location $WorkDir
try {
    $ScriptsDir = (Resolve-Path "..\scripts\analysis").Path

    Write-Host "========== Starting Paper Analysis for System: $System ========" -ForegroundColor Green
    Write-Host ">> Using Python Interpreter: $PY" -ForegroundColor Green

    Write-Host "[0/8] Generating analysis index groups (0_make_index.sh) ..."
    & bash "$ScriptsDir/0_make_index.sh"

    Write-Host "[1/8] Backbone RMSD / RMSF (1_rmsd_rmsf.sh) ..."
    & bash "$ScriptsDir/1_rmsd_rmsf.sh"

    Write-Host "[2/8] Peptide around AChE COM RDF (2_rdf.sh) ..."
    & bash "$ScriptsDir/2_rdf.sh"

    Write-Host "[3/8] Complex SASA (3_sasa.sh) ..."
    & bash "$ScriptsDir/3_sasa.sh"

    Write-Host "[4/8] Peptide Secondary Structure Evolution (4_secondary_structure.sh) ..."
    & bash "$ScriptsDir/4_secondary_structure.sh"

    Write-Host "[5/8] Inter/Intra Hydrogen Bonds (5_hbond.sh) ..."
    & bash "$ScriptsDir/5_hbond.sh"

    Write-Host "[6/8] Native/Non-native Intermolecular Contacts (contacts.py) ..."
    & $PY "$ScriptsDir\contacts.py" -t md.tpr -f md.xtc

    Write-Host "[7/8] Water-mediated Bridging Interactions (bridging_waters.py) ..."
    & $PY "$ScriptsDir\bridging_waters.py" -t md.tpr -f md.xtc

    Write-Host "[8/8] Generating publication SVG / PNG / PDF figures (plot_all.py) ..."
    & $PY "$ScriptsDir\plot_all.py" --dir . --out ./figures

    Write-Host "========== System $System Analysis Finished Successfully! ==========" -ForegroundColor Green
    Write-Host "Generated Figures in ./figures/ :"
    Write-Host "  - fig1_rmsd_rmsf.{svg,png,pdf}            -> Figure 1 (RMSD / RMSF)"
    Write-Host "  - fig2_rdf.{svg,png,pdf}                  -> Figure 2 (RDF)"
    Write-Host "  - fig3_sasa.{svg,png,pdf}                 -> Figure 3 (SASA)"
    Write-Host "  - fig4_secondary_structure.{svg,png,pdf}  -> Figure 4 (DSSP Fraction)"
    Write-Host "  - fig5_contacts.{svg,png,pdf}             -> Figure 5 (Contacts)"
    Write-Host "  - fig6_bridging_waters.{svg,png,pdf}      -> Figure 6 (Bridging Waters)"
    Write-Host "  - fig_hbonds.{svg,png,pdf}                -> Paper 3.3 (H-Bonds)"
    Write-Host "  - fig0_summary_all.{svg,png,pdf}          -> 2x3 Combined Master Plot"
    Write-Host "  - summary_metrics.csv and wide.csv        -> Summary Tables"
} finally {
    Pop-Location
}
