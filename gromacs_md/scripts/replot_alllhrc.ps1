# Replot ONE finished v1.0 system. Does not start mdrun. Does not delete md_*.
# Default: OnlyPlot if analysis products already exist (safe after a plot crash).
#   .\replot_alllhrc.ps1
#   .\replot_alllhrc.ps1 -OnlyPlot
#   .\replot_alllhrc.ps1 -Full
param(
    [switch]$OnlyPlot,
    [switch]$Full
)
try { chcp 65001 | Out-Null } catch {}
$ErrorActionPreference = "Continue"
Write-Host "========== REPLOT alllhrc ==========" -ForegroundColor Cyan
& "$PSScriptRoot\replot_common.ps1" -System alllhrc -OnlyPlot:$OnlyPlot -Full:$Full
exit $LASTEXITCODE
