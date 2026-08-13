# ============================================================
# AChE - A-beta Complex Trajectory Analysis & Figure Generation (PowerShell)
#
# Usage:  .\run_analysis.ps1 -System alllhrc
#         .\run_analysis.ps1 -System alllhrc -Testing
#         .\run_analysis.ps1 -System alllhrc -OnlyPlot
#         .\run_analysis.ps1 -System alllhrc -OnlyDsspPeptide   # 只修 DSSP+肽三段, 不重跑其它分析
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System,

    [switch]$Testing,
    [switch]$OnlyPlot,
    [switch]$OnlyDsspPeptide
)

[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::InputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

$GMX = "gmx"
if (Get-Command "gmx.exe" -ErrorAction SilentlyContinue) {
    $GMX = "gmx.exe"
}

$PY = "python"
if (Get-Command "python.exe" -ErrorAction SilentlyContinue) {
    $PY = "python.exe"
}

$WorkDir = Join-Path (Split-Path $PSScriptRoot -Parent) "md_$System"
if (-not (Test-Path $WorkDir)) {
    $WorkDir = "..\md_$System"
}
try { chcp 65001 | Out-Null } catch {}

if (-not (Test-Path "$WorkDir\md.xtc") -and -not (Test-Path "$WorkDir\md_0_1.xtc")) {
    Write-Host "ERROR: no md.xtc / md_0_1.xtc in $WorkDir" -ForegroundColor Red
    Write-Host "This script is ANALYSIS only. Do NOT run run_all_four_100ns_formal.ps1 (it deletes md_*)." -ForegroundColor Yellow
    exit 1
}

$TprFile = if (Test-Path "$WorkDir\md.tpr") { "md.tpr" } else { "md_0_1.tpr" }
$XtcFile = if (Test-Path "$WorkDir\md.xtc") { "md.xtc" } else { "md_0_1.xtc" }

if ($Testing -or $env:TESTING -eq "1") {
    Write-Host ">> [TESTING MODE] 5000-step short trajectory parameters" -ForegroundColor Yellow
    $env:TESTING = "1"
} else {
    Write-Host ">> [PRODUCTION MODE] 100 ns production trajectory analysis" -ForegroundColor Cyan
    $env:TESTING = "0"
}

Push-Location $WorkDir
try {
    Write-Host "========== Starting Paper Analysis for System: $System ========" -ForegroundColor Green
    Write-Host ">> Using Python Interpreter: $PY" -ForegroundColor Green

    if ($OnlyDsspPeptide) {
        Write-Host ">> [OnlyDsspPeptide] 只重跑 DSSP + 肽三段 + 出图, 其它分析全部跳过" -ForegroundColor Yellow
        if (-not (Test-Path "index.ndx")) {
            Write-Host "[0] index.ndx 缺失, 补建 ..."
            & bash "../scripts/analysis/0_make_index.sh" 2>&1 | ForEach-Object { "$_" }
        }
        Write-Host "[DSSP] 4_secondary_structure.sh ..."
        & bash "../scripts/analysis/4_secondary_structure.sh" 2>&1 | ForEach-Object { "$_" }
        Write-Host "[Peptide phases] analyze_peptide_phases.py ..."
        & $PY "..\scripts\analysis\analyze_peptide_phases.py" -d . 2>&1 | ForEach-Object { "$_" }
    } elseif (-not $OnlyPlot) {
        Write-Host "[0/8] Generating analysis index groups (0_make_index.sh) ..."
        & bash "../scripts/analysis/0_make_index.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[1/8] Backbone RMSD / RMSF (1_rmsd_rmsf.sh) ..."
        & bash "../scripts/analysis/1_rmsd_rmsf.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[2/8] Peptide around AChE COM RDF (2_rdf.sh) ..."
        & bash "../scripts/analysis/2_rdf.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[3/8] Complex SASA (3_sasa.sh) ..."
        & bash "../scripts/analysis/3_sasa.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[4/9] Complex DSSP (4_secondary_structure.sh) ..."
        & bash "../scripts/analysis/4_secondary_structure.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[5/9] Radius of gyration Rg (6_rg.sh) ..."
        & bash "../scripts/analysis/6_rg.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[6/9] Inter/Intra Hydrogen Bonds (5_hbond.sh) ..."
        & bash "../scripts/analysis/5_hbond.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[7/9] Native/Non-native Intermolecular Contacts (contacts.py) ..."
        & $PY "..\scripts\analysis\contacts.py" -t $TprFile -f $XtcFile 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[7/8] Water-mediated Bridging Interactions (bridging_waters.py) ..."
        & $PY "..\scripts\analysis\bridging_waters.py" -t $TprFile -f $XtcFile 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    } else {
        Write-Host ">> [OnlyPlot Mode] Skipping GROMACS analysis steps [0-7], directly generating figures via plot_all.py ..." -ForegroundColor Yellow
    }

    Write-Host "[9/9] Generating publication SVG / PNG / PDF figures (plot_all.py) ..."
    & $PY "..\scripts\analysis\plot_all.py" --dir . --out ./figures 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    Write-Host "========== System $System Analysis Finished Successfully! ==========" -ForegroundColor Green
    Write-Host "Generated Figures in ./figures/ :"
    Write-Host "  - fig1_rmsd_rmsf.{svg,png,pdf}            => Complex/AChE RMSD + AChE RMSF"
    Write-Host "  - fig_peptide_rmsd_rmsf.{svg,png,pdf}     => Peptide RMSD + RMSF (standalone)"
    Write-Host "  - fig2_rdf.{svg,png,pdf}                  => Figure 2 (RDF)"
    Write-Host "  - fig3_sasa.{svg,png,pdf}                 => Figure 3 (SASA)"
    Write-Host "  - fig4_secondary_structure.{svg,png,pdf}  => DSSP % lines + occupancy bars + heatmap + peptide"
    Write-Host "  - fig5_contacts.{svg,png,pdf}             => Figure 5 (Contacts)"
    Write-Host "  - fig6_bridging_waters.{svg,png,pdf}      => Figure 6 (Bridging Waters)"
    Write-Host "  - fig_hbonds.{svg,png,pdf}                => H-Bonds (standalone)"
    Write-Host "  - fig0_summary_all.{svg,png,pdf}          => 2x4: RMSD RMSF RDF SASA | DSSP% Rg H-bonds DSSP-lines"
    Write-Host "  - summary_metrics.csv & wide.csv          => Summary Tables"
} finally {
    Pop-Location
}
