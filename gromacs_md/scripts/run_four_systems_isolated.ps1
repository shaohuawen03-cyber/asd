param([switch]$Testing)
$ErrorActionPreference = "Stop"
$SourceRoot = Split-Path $PSScriptRoot -Parent
$stamp = Get-Date -Format "yyyyMMdd-HHmmss"
$RunRoot = Join-Path $SourceRoot ("rerun_four_systems-" + $stamp)
$RunScripts = Join-Path $RunRoot "scripts"
New-Item -ItemType Directory -Path $RunRoot -Force | Out-Null
Write-Host "NEW isolated workspace: $RunRoot" -ForegroundColor Green
Copy-Item (Join-Path $SourceRoot "input") $RunRoot -Recurse
Copy-Item (Join-Path $SourceRoot "mdp") $RunRoot -Recurse
Copy-Item (Join-Path $SourceRoot "scripts") $RunRoot -Recurse
function Run-Step([string]$label, [scriptblock]$action) {
  Write-Host "=== $label ===" -ForegroundColor Cyan
  & $action
  if ($LASTEXITCODE -ne 0) { throw "$label failed with exit code $LASTEXITCODE" }
}
Push-Location $RunScripts
try {
  foreach ($s in @("alllhrc","fllhttr","ylsllqr")) {
    if ($Testing) { Run-Step "MD $s testing" { & .\run_all.ps1 -System $s -Testing } }
    else { Run-Step "MD $s 100 ns" { & .\run_all.ps1 -System $s } }
    Run-Step "Analysis $s" { & .\run_analysis.ps1 -System $s }
  }
  if ($Testing) { Run-Step "MD and analysis ache apo testing" { & .\run_apo_ache_100ns.ps1 -Testing } }
  else { Run-Step "MD and analysis ache apo 100 ns" { & .\run_apo_ache_100ns.ps1 } }
  foreach ($s in @("ache","alllhrc","fllhttr","ylsllqr")) {
    $p = Join-Path $RunRoot "md_$s\figures\summary_metrics.csv"
    if (-not (Test-Path $p)) { throw "Missing final summary: $p" }
  }
  Write-Host "FINAL VALIDATION PASSED: all four summary_metrics.csv exist" -ForegroundColor Green
} finally { Pop-Location }
Write-Host "DONE: $RunRoot" -ForegroundColor Green
