# Replot an existing isolated rerun only after validating all four analysis inputs.
# Usage: .\replot_rerun_validated.ps1 -RunRoot ..\rerun_four_systems-YYYYMMDD-HHMMSS
param([Parameter(Mandatory=$true)][string]$RunRoot)
$ErrorActionPreference = "Stop"
$RunRoot = (Resolve-Path $RunRoot).Path
$systems = @("ache","alllhrc","fllhttr","ylsllqr")
foreach ($s in $systems) {
    $d = Join-Path $RunRoot "md_$s"
    foreach ($f in @("md.tpr","md.xtc")) {
        if (-not (Test-Path (Join-Path $d $f))) { throw "Missing $d\$f; do not plot." }
    }
    $summary = Join-Path $d "figures\summary_metrics.csv"
    if (-not (Test-Path $summary)) { throw "Missing $summary; analysis did not finish." }
    $rows = @(Import-Csv $summary)
    if ($rows.Count -eq 0) { throw "Empty $summary; do not plot." }
    Write-Host "OK ${s}: trajectory and summary present ($($rows.Count) rows)" -ForegroundColor Green
}
$scriptDir = Join-Path $RunRoot "scripts"
Push-Location $scriptDir
try {
    & .\run_unified_replot.ps1
    if ($LASTEXITCODE -ne 0) { throw "run_unified_replot.ps1 failed with exit code $LASTEXITCODE" }
} finally { Pop-Location }
Write-Host "Validated plots are in $RunRoot\compare_ache_vs_*" -ForegroundColor Green
