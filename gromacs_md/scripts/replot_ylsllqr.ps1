# Replot ONE finished v1.0 system. Does not start mdrun. Does not delete md_*.
# Default: OnlyPlot if analysis products already exist (safe after a plot crash).
#   .\replot_ylsllqr.ps1
#   .\replot_ylsllqr.ps1 -OnlyPlot
#   .\replot_ylsllqr.ps1 -Full
param(
    [switch]$OnlyPlot,
    [switch]$Full
)
try { chcp 65001 | Out-Null } catch {}
$ErrorActionPreference = "Continue"
Write-Host "========== REPLOT ylsllqr ==========" -ForegroundColor Cyan
& "$PSScriptRoot\replot_common.ps1" -System ylsllqr -OnlyPlot:$OnlyPlot -Full:$Full
exit $LASTEXITCODE
