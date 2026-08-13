# Re-analyze + replot ONE finished v1.0 system. Does not start mdrun.
# DSSP = complex (v2.5). fig0 = no peptide RMSD/RMSF; adds complex Rg.
param()
$ErrorActionPreference = "Continue"
try { chcp 65001 | Out-Null } catch {}
Write-Host "========== REPLOT alllhrc (v1.0 traj + current analysis) ==========" -ForegroundColor Cyan
& "$PSScriptRoot\run_analysis.ps1" -System alllhrc
exit $LASTEXITCODE
