# Re-analyze / replot finished 100 ns systems with the CURRENT analysis code.
# Does NOT delete md_* directories. Does NOT start mdrun.
#
# Usage:
#   .\rerun_analysis_four.ps1
#   .\rerun_analysis_four.ps1 -Systems alllhrc,fllhttr
#   .\rerun_analysis_four.ps1 -OnlyPlot
#   .\rerun_analysis_four.ps1 -Full
param(
    [string[]]$Systems = @("alllhrc", "fllhttr", "ylsllqr", "ache"),
    [switch]$OnlyPlot,
    [switch]$Full
)

try { chcp 65001 | Out-Null } catch {}
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

Write-Host "============================================================" -ForegroundColor Green
Write-Host " Replot finished v1.0 systems (PBC md_fit + complex DSSP + plots)" -ForegroundColor Green
Write-Host " Will NOT delete trajectories. Will NOT start mdrun." -ForegroundColor Yellow
Write-Host "============================================================" -ForegroundColor Green

$Results = @()
$Fail = 0
foreach ($Sys in $Systems) {
    Write-Host ""
    & "$PSScriptRoot\replot_common.ps1" -System $Sys -OnlyPlot:$OnlyPlot -Full:$Full
    $code = $LASTEXITCODE
    if ($code -eq 2) {
        $Results += [PSCustomObject]@{ System = $Sys; Status = "SKIP_MDRUN_STILL_RUNNING" }
    } elseif ($code -eq 1) {
        $Results += [PSCustomObject]@{ System = $Sys; Status = "SKIP_NO_XTC_OR_DIR" }
    } elseif ($code -ne 0 -and $null -ne $code) {
        Write-Host "[FAIL] $Sys exit $code" -ForegroundColor Red
        $Results += [PSCustomObject]@{ System = $Sys; Status = "FAIL_$code" }
        $Fail = 1
    } else {
        $Results += [PSCustomObject]@{ System = $Sys; Status = "OK" }
    }
}

Write-Host ""
Write-Host "============================================================" -ForegroundColor Green
$Results | Format-Table -AutoSize
Write-Host "Do NOT run these again (they DELETE md_* and restart MD):"
Write-Host "  run_all_four_100ns_formal.ps1 / run_100ns_formal.ps1 / run_all.sh"
Write-Host "  run_split_md_workflow.ps1 / run_all_four_user_workflow.ps1"
Write-Host "  clean_test_results.ps1"
Write-Host "============================================================" -ForegroundColor Green
if ($Fail -ne 0) { exit 1 }
exit 0
