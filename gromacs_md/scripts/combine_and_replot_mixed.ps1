# Build a new comparison workspace without overwriting any prior results.
# Uses first/original alllhrc plus the final successful rerun apo/fllhttr/ylsllqr.
param(
    [Parameter(Mandatory=$true)][string]$FinalRunRoot,
    [string]$OldAlllhrc = "..\md_alllhrc"
)
$ErrorActionPreference = "Stop"
$Here = $PSScriptRoot
$Root = Split-Path $Here -Parent
$FinalRunRoot = (Resolve-Path $FinalRunRoot).Path
$OldAlllhrc = (Resolve-Path (Join-Path $Here $OldAlllhrc)).Path
$stamp = Get-Date -Format "yyyyMMdd-HHmmss"
$Out = Join-Path $Root ("combined_compare_mixed-" + $stamp)
New-Item -ItemType Directory -Force $Out | Out-Null
Write-Host "New comparison workspace: $Out" -ForegroundColor Green

# Copy workflow code and inputs needed by the plotting scripts.
Copy-Item (Join-Path $Root "scripts") (Join-Path $Out "scripts") -Recurse
Copy-Item (Join-Path $Root "mdp") (Join-Path $Out "mdp") -Recurse
Copy-Item (Join-Path $Root "input") (Join-Path $Out "input") -Recurse

# First/original alllhrc result; final rerun result for the other three systems.
Copy-Item $OldAlllhrc (Join-Path $Out "md_alllhrc") -Recurse
foreach ($s in @("ache","fllhttr","ylsllqr")) {
    $src = Join-Path $FinalRunRoot ("md_" + $s)
    if (-not (Test-Path $src)) { throw "Missing final rerun directory: $src" }
    Copy-Item $src (Join-Path $Out ("md_" + $s)) -Recurse
}

# Validate before plotting. No blank comparison is allowed.
foreach ($s in @("ache","alllhrc","fllhttr","ylsllqr")) {
    $d = Join-Path $Out ("md_" + $s)
    foreach ($f in @("md.tpr","md.xtc")) {
        if (-not (Test-Path (Join-Path $d $f))) { throw "Missing $d\$f" }
    }
    $summary = Join-Path $d "figures\summary_metrics.csv"
    if (-not (Test-Path $summary)) { throw "Missing $summary" }
    if (@(Import-Csv $summary).Count -eq 0) { throw "Empty $summary" }
}

Push-Location (Join-Path $Out "scripts")
try {
    & .\run_unified_replot.ps1
    if ($LASTEXITCODE -ne 0) { throw "run_unified_replot.ps1 failed: $LASTEXITCODE" }
} finally { Pop-Location }

foreach ($s in @("alllhrc","fllhttr","ylsllqr")) {
    $p = Join-Path $Out ("compare_ache_vs_" + $s + "\compare_summary.csv")
    if (-not (Test-Path $p)) { throw "Missing comparison output: $p" }
}
Write-Host "SUCCESS: comparison outputs are in $Out\compare_ache_vs_*" -ForegroundColor Green
Write-Host "Original and all rerun directories were not modified." -ForegroundColor Green
