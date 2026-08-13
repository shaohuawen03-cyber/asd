# Re-analyze finished 100 ns systems with the CURRENT analysis code.
# Does NOT delete md_* directories. Does NOT start mdrun.
#
# Usage:
#   .\rerun_analysis_four.ps1
#   .\rerun_analysis_four.ps1 -Systems alllhrc,fllhttr
#
param(
    [string[]]$Systems = @("alllhrc", "fllhttr", "ylsllqr", "ache")
)

try { chcp 65001 | Out-Null } catch {}
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

Write-Host "============================================================" -ForegroundColor Green
Write-Host " Re-run ANALYSIS only (PBC md_fit + DSSP + plots)" -ForegroundColor Green
Write-Host " Will NOT delete trajectories. Will NOT start mdrun." -ForegroundColor Yellow
Write-Host "============================================================" -ForegroundColor Green

$Results = @()
foreach ($Sys in $Systems) {
    $Work = "..\md_$Sys"
    $Xtc = $null
    foreach ($c in @("$Work\md_0_1.xtc", "$Work\md.xtc")) {
        if (Test-Path $c) { $Xtc = $c; break }
    }
    if (-not $Xtc) {
        Write-Host "[SKIP] $Sys : no md_0_1.xtc / md.xtc (MD not finished or not started)" -ForegroundColor Yellow
        $Results += [PSCustomObject]@{ System = $Sys; Status = "SKIP_NO_XTC" }
        continue
    }

    # gro is written when mdrun finishes; if missing, production is probably still running
    $GroDone = (Test-Path "$Work\md_0_1.gro") -or (Test-Path "$Work\md.gro")
    if (-not $GroDone) {
        Write-Host "[SKIP] $Sys : xtc exists but no md_0_1.gro -- mdrun likely still running" -ForegroundColor Yellow
        $Results += [PSCustomObject]@{ System = $Sys; Status = "SKIP_MDRUN_STILL_RUNNING" }
        continue
    }

    Write-Host ""
    Write-Host "========== ANALYZE $Sys ==========" -ForegroundColor Cyan
    & .\run_analysis.ps1 -System $Sys
    $code = $LASTEXITCODE
    if ($code -ne 0 -and $null -ne $code) {
        Write-Host "[FAIL] $Sys analysis exit $code" -ForegroundColor Red
        $Results += [PSCustomObject]@{ System = $Sys; Status = "FAIL_$code" }
    } else {
        Write-Host "[OK] $Sys -> md_$Sys\figures\" -ForegroundColor Green
        $Results += [PSCustomObject]@{ System = $Sys; Status = "OK" }
    }
}

Write-Host ""
Write-Host "============================================================" -ForegroundColor Green
$Results | Format-Table -AutoSize
Write-Host "Do NOT run run_split_md_workflow.ps1 or run_all_four_user_workflow.ps1"
Write-Host "again: those scripts DELETE md_* and restart MD from zero."
Write-Host "============================================================" -ForegroundColor Green
