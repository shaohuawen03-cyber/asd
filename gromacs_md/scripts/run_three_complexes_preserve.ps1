# Run the three AChE-peptide complexes in a NEW timestamped workspace.
# The original md_* and compare_* directories are never touched.
param([switch]$Testing)
$ErrorActionPreference = "Stop"
$SourceRoot = Split-Path $PSScriptRoot -Parent
$stamp = Get-Date -Format "yyyyMMdd-HHmmss"
$RunRoot = Join-Path $SourceRoot ("rerun_complexes-" + $stamp)
$RunScripts = Join-Path $RunRoot "scripts"
New-Item -ItemType Directory -Path $RunRoot -Force | Out-Null

Write-Host "New isolated workspace: $RunRoot" -ForegroundColor Green
# Copy only workflow inputs and code; no previous md_* results or trajectories are copied.
Copy-Item (Join-Path $SourceRoot "input") -Destination $RunRoot -Recurse
Copy-Item (Join-Path $SourceRoot "mdp") -Destination $RunRoot -Recurse
Copy-Item (Join-Path $SourceRoot "scripts") -Destination $RunRoot -Recurse

$systems = @("alllhrc", "fllhttr", "ylsllqr")
Push-Location $RunScripts
try {
    foreach ($s in $systems) {
        Write-Host "=== Running $s in isolated workspace ===" -ForegroundColor Cyan
        if ($Testing) { & .\run_all.ps1 -System $s -Testing }
        else { & .\run_all.ps1 -System $s }
        if ($LASTEXITCODE -ne 0) { throw "run_all.ps1 failed for $s with exit code $LASTEXITCODE" }
    }
    Write-Host "=== Generating final comparison figures ===" -ForegroundColor Cyan
    & .\run_unified_replot.ps1
    if ($LASTEXITCODE -ne 0) { throw "run_unified_replot.ps1 failed with exit code $LASTEXITCODE" }
} finally { Pop-Location }
Write-Host "DONE. New results: $RunRoot\md_* and $RunRoot\compare_ache_vs_*" -ForegroundColor Green
Write-Host "Original results were not modified." -ForegroundColor Green
