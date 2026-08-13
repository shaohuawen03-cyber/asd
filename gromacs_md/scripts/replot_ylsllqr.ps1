# Re-analyze + replot ONE finished v1.0 system. Does not start mdrun.
param()
$ErrorActionPreference = "Continue"
try { chcp 65001 | Out-Null } catch {}
Write-Host "========== REPLOT ylsllqr (v1.0 traj + current analysis) ==========" -ForegroundColor Cyan
& "$PSScriptRoot\run_analysis.ps1" -System ylsllqr
exit $LASTEXITCODE
