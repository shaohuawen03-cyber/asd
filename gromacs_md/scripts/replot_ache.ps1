# Replot AChE monomer. Does not start mdrun. Does not delete md_*.
# Default: OnlyPlot if analysis products already exist (safe after a plot crash).
#   .\replot_ache.ps1
#   .\replot_ache.ps1 -OnlyPlot
#   .\replot_ache.ps1 -Full
param(
    [switch]$OnlyPlot,
    [switch]$Full
)
try { chcp 65001 | Out-Null } catch {}
$ErrorActionPreference = "Continue"
Write-Host "========== REPLOT ache monomer ==========" -ForegroundColor Cyan
& "$PSScriptRoot\replot_common.ps1" -System ache -OnlyPlot:$OnlyPlot -Full:$Full
exit $LASTEXITCODE
