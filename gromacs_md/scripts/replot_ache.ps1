# Re-analyze + replot AChE monomer. Does not start mdrun.
param()
$ErrorActionPreference = "Continue"
try { chcp 65001 | Out-Null } catch {}
Write-Host "========== REPLOT ache monomer (v1.0 traj + current analysis) ==========" -ForegroundColor Cyan
& "$PSScriptRoot\run_analysis.ps1" -System ache
exit $LASTEXITCODE
