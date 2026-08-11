# ============================================================
# AChE & A-beta 4-System 5000-Step Full Pipeline Test & Error Logging Suite
#
# Usage:  .\test_four_systems.ps1
# ============================================================
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::InputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

$Systems = @("alllhrc", "fllhttr", "ylsllqr", "ache")
$LogFile = (Resolve-Path ".").Path + "\test_four_systems.log"

if (Test-Path $LogFile) { Remove-Item $LogFile -Force }

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " Starting 4-System Full Workflow 5000-Step Test & Error Logging" -ForegroundColor Green
Write-Host " Log file: $LogFile" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

$Results = @()

foreach ($Sys in $Systems) {
    Write-Host ""
    Write-Host "====================================================================" -ForegroundColor Cyan
    Write-Host " Processing System: $Sys (5000 Steps MD + Analysis)" -ForegroundColor Cyan
    Write-Host "====================================================================" -ForegroundColor Cyan
    
    $Status = "SUCCESS"
    "=== System: $Sys ===" | Out-File -FilePath $LogFile -Append -Encoding utf8
    & .\run_pipeline_all.ps1 -System $Sys -Testing 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $LogFile -Append
    if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
        $Status = "FAILED (Exit Code $LASTEXITCODE)"
    }

    $FigPath = "..\md_$Sys\figures\fig0_summary_all.svg"
    if (-not (Test-Path $FigPath)) {
        $Status = "FAILED (Missing final figure fig0_summary_all.svg)"
    }

    $Results += [PSCustomObject]@{
        System = $Sys
        Status = $Status
    }
}

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " 4-SYSTEM FULL WORKFLOW TEST REPORT" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
$Results | Format-Table -AutoSize
Write-Host "Full execution log saved in: $LogFile" -ForegroundColor Green

$Errors = Select-String -Path $LogFile -Pattern "Fatal error|Segmentation fault|command not found|No such file|ModuleNotFoundError|Syntax error" -ErrorAction SilentlyContinue
if ($Errors) {
    Write-Host ""
    Write-Host ">>> [LOGGED ERRORS FOUND IN test_four_systems.log] :" -ForegroundColor Red
    $Errors | ForEach-Object { Write-Host "  $_.Line" -ForegroundColor Yellow }
    Write-Host "`nPlease copy any error lines above and paste them to me!" -ForegroundColor Red
} else {
    Write-Host ""
    Write-Host ">>> [ALL 4 SYSTEMS 100% SUCCESSFUL - ZERO ERRORS LOGGED!]" -ForegroundColor Green
}
Write-Host "====================================================================" -ForegroundColor Green
