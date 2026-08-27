# Re-run the three AChE-peptide complexes without overwriting previous results.
# Run from gromacs_md\scripts. Existing md_* directories are moved to a timestamped archive.
param([switch]$Testing)
$ErrorActionPreference = "Stop"
$Here = $PSScriptRoot
$Root = Split-Path $Here -Parent
$stamp = Get-Date -Format "yyyyMMdd-HHmmss"
$Archive = Join-Path $Root ("rerun_archive_complexes-" + $stamp)
New-Item -ItemType Directory -Path $Archive -Force | Out-Null
$systems = @("alllhrc", "fllhttr", "ylsllqr")
Write-Host "Archive: $Archive" -ForegroundColor Yellow

foreach ($s in $systems) {
    $work = Join-Path $Root ("md_" + $s)
    if (Test-Path $work) {
        $dest = Join-Path $Archive ("md_" + $s)
        Write-Host "Backing up $work -> $dest" -ForegroundColor Yellow
        Move-Item -LiteralPath $work -Destination $dest
    }
    $compare = Join-Path $Root ("compare_ache_vs_" + $s)
    if (Test-Path $compare) {
        $destCompare = Join-Path $Archive ("compare_ache_vs_" + $s)
        Write-Host "Backing up $compare -> $destCompare" -ForegroundColor Yellow
        Move-Item -LiteralPath $compare -Destination $destCompare
    }
}

Push-Location $Here
try {
    foreach ($s in $systems) {
        Write-Host "=== Running $s ===" -ForegroundColor Cyan
        if ($Testing) {
            & .\run_all.ps1 -System $s -Testing
        } else {
            & .\run_all.ps1 -System $s
        }
        if ($LASTEXITCODE -ne 0) { throw "run_all.ps1 failed for $s with exit code $LASTEXITCODE" }
    }
    Write-Host "=== Replotting three complexes and comparisons ===" -ForegroundColor Cyan
    & .\run_unified_replot.ps1
    if ($LASTEXITCODE -ne 0) { throw "run_unified_replot.ps1 failed with exit code $LASTEXITCODE" }
} finally { Pop-Location }
Write-Host "DONE. New results are in md_alllhrc, md_fllhttr, md_ylsllqr and compare_*; old results are preserved in $Archive" -ForegroundColor Green
