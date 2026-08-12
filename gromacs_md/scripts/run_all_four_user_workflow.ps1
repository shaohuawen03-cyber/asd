# ============================================================
# 4-System 100 ns MD & Analysis Master Sequential Suite (User Protocol)
#
# Automatically runs 100 ns simulation & plotting for all 4 systems sequentially:
#   1. Clean up all previous test directories and log files
#   2. System 1 (alllhrc): 100 ns MD (md_0_1) -> full analysis -> figures ready for inspection
#   3. System 2 (fllhttr): 100 ns MD (md_0_1) -> full analysis -> figures ready for inspection
#   4. System 3 (ylsllqr): 100 ns MD (md_0_1) -> full analysis -> figures ready for inspection
#   5. System 4 (ache monomer): 100 ns MD (md_0_1) -> full analysis -> figures ready for inspection
#   6. Comparative plotting: ache vs each complex
#
# Usage:  .\run_all_four_user_workflow.ps1
# ============================================================
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::InputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

# Force clear any residual testing mode environment variables
$env:TESTING = "0"
$env:ONLY_MD = "0"
Remove-Item env:TESTING -ErrorAction SilentlyContinue
Remove-Item env:ONLY_MD -ErrorAction SilentlyContinue

$Systems = @("alllhrc", "fllhttr", "ylsllqr", "ache")
$LogFile = (Resolve-Path ".").Path + "\run_all_four_user_workflow.log"

if (Test-Path $LogFile) { Remove-Item $LogFile -Force }

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " [Step 0] Cleaning up previous test results to ensure 100% clean start..." -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
.\clean_test_results.ps1

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " Starting Sequential 100 ns MD & Analysis for 4 Systems (User Workflow)" -ForegroundColor Green
Write-Host " Systems in queue: $($Systems -join ', ')" -ForegroundColor Green
Write-Host " Log file: $LogFile" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

$Results = @()
$StartTime = Get-Date

foreach ($Sys in $Systems) {
    $SysStart = Get-Date
    Write-Host ""
    Write-Host "====================================================================" -ForegroundColor Cyan
    Write-Host " >>> [STARTING 100 NS SIMULATION & ANALYSIS] System: $Sys" -ForegroundColor Cyan
    Write-Host "     Timestamp: $($SysStart.ToString('yyyy-MM-dd HH:mm:ss'))" -ForegroundColor Cyan
    Write-Host "====================================================================" -ForegroundColor Cyan

    $Status = "SUCCESS"
    "=== [100NS START] System: $Sys at $($SysStart.ToString('yyyy-MM-dd HH:mm:ss')) ===" | Out-File -FilePath $LogFile -Append -Encoding utf8

    & .\run_md_user_workflow.ps1 -System $Sys 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $LogFile -Append
    if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
        $Status = "FAILED_MD (Exit Code $LASTEXITCODE)"
    } else {
        & .\run_analysis.ps1 -System $Sys 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $LogFile -Append
        if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
            $Status = "FAILED_ANALYSIS (Exit Code $LASTEXITCODE)"
        }
    }

    $FigPath = "..\md_$Sys\figures\fig0_summary_all.svg"
    if (-not (Test-Path $FigPath)) {
        $Status = "FAILED (Missing final figure fig0_summary_all.svg)"
    }

    $SysEnd = Get-Date
    $Duration = ($SysEnd - $SysStart).ToString("dd\.hh\:mm\:ss")
    
    $Results += [PSCustomObject]@{
        System   = $Sys
        Status   = $Status
        Duration = $Duration
    }

    Write-Host ""
    Write-Host "====================================================================" -ForegroundColor Green
    Write-Host " [READY FOR YOUR INSPECTION] System $Sys 100 ns MD and Figures Completed!" -ForegroundColor Green
    Write-Host " -> Open directory: md_$Sys/figures/ to inspect SVG/PNG/PDF plots & CSV tables!" -ForegroundColor Green
    Write-Host " -> Moving directly to next system in queue without delay..." -ForegroundColor Green
    Write-Host "====================================================================" -ForegroundColor Green
}

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Cyan
Write-Host " >>> [BONUS STEP] Generating AChE Monomer vs Complex Comparative Figures..." -ForegroundColor Cyan
Write-Host "====================================================================" -ForegroundColor Cyan
foreach ($CompSys in @("alllhrc", "fllhttr", "ylsllqr")) {
    if ((Test-Path "..\md_ache\md_0_1.xtc") -and (Test-Path "..\md_$CompSys\md_0_1.xtc")) {
        Write-Host "  -> Comparing ache Monomer vs $CompSys Complex..." -ForegroundColor Yellow
        & .\plot_compare.ps1 -Protein "ache" -Complex $CompSys -OutDir "..\compare_ache_vs_$CompSys" 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $LogFile -Append
    }
}

$TotalTime = (Get-Date - $StartTime).ToString("dd\.hh\:mm\:ss")

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " 4-SYSTEM 100 NS MD & ANALYSIS MASTER REPORT (User Workflow)" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
$Results | Format-Table -AutoSize
Write-Host "Total Pipeline Duration: $TotalTime" -ForegroundColor Green
Write-Host "Full execution log saved in: $LogFile" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
