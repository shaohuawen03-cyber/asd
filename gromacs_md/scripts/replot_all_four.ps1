# One-click replot of all finished v1.0 systems.
# Does NOT start mdrun. Does NOT delete md_*.
#
#   .\replot_all_four.ps1
#   .\replot_all_four.ps1 -OnlyPlot
#   .\replot_all_four.ps1 -Full
param(
    [switch]$OnlyPlot,
    [switch]$Full
)
try { chcp 65001 | Out-Null } catch {}
& "$PSScriptRoot\rerun_analysis_four.ps1" -OnlyPlot:$OnlyPlot -Full:$Full
exit $LASTEXITCODE
