# Replot ONE finished v1.0 system. Does not start mdrun. Does not delete md_*.
# Default: OnlyPlot if analysis products already exist (safe after a plot crash).
#   .\replot_fllhttr.ps1
#   .\replot_fllhttr.ps1 -OnlyPlot
#   .\replot_fllhttr.ps1 -Full
param(
    [switch]$OnlyPlot,
    [switch]$Full
)
try { chcp 65001 | Out-Null } catch {}
$ErrorActionPreference = "Continue"
Write-Host "========== REPLOT fllhttr ==========" -ForegroundColor Cyan
& "$PSScriptRoot\replot_common.ps1" -System fllhttr -OnlyPlot:$OnlyPlot -Full:$Full
exit $LASTEXITCODE
