# One-click replot of all finished v1.0 systems.
# Same as rerun_analysis_four.ps1. Does NOT start mdrun. Does NOT delete md_*.
& "$PSScriptRoot\rerun_analysis_four.ps1" @args
exit $LASTEXITCODE
