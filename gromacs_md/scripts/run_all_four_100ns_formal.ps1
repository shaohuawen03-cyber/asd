# ============================================================
# 4-System 100 ns Formal Production MD & Analysis Sequential Suite (PowerShell)
#
# Automatically runs 100 ns formal simulation & plotting for all 4 systems sequentially:
#   1. Clean up all previous test directories and log files
#   2. System 1 (alllhrc): 100 ns MD -> full analysis -> figures ready for inspection
#   3. System 2 (fllhttr): 100 ns MD -> full analysis -> figures ready for inspection
#   4. System 3 (ylsllqr): 100 ns MD -> full analysis -> figures ready for inspection
#   5. System 4 (ache monomer): 100 ns MD -> full analysis -> figures ready for inspection
#   6. Comparative plotting: ache vs each complex
#
# Usage:  .\run_all_four_100ns_formal.ps1
# ============================================================
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::InputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

# 强制清除之前任何测试会话残留的环境变量，1000% 锁定为 100 ns 正式产物生产模式
$env:TESTING = "0"
$env:ONLY_MD = "0"
$env:SIM_MODE = "100ns"
Remove-Item env:TESTING -ErrorAction SilentlyContinue
Remove-Item env:ONLY_MD -ErrorAction SilentlyContinue

$Systems = @("alllhrc", "fllhttr", "ylsllqr", "ache")
$LogFile = (Resolve-Path ".").Path + "\run_all_four_100ns_formal.log"

if (Test-Path $LogFile) { Remove-Item $LogFile -Force }

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " [Step 0] Cleaning up previous test results to ensure 100% clean formal start..." -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
.\clean_test_results.ps1

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " Starting Sequential 100 ns Formal Production MD & Analysis for 4 Systems" -ForegroundColor Green
Write-Host " Systems in queue: $($Systems -join ', ')" -ForegroundColor Green
Write-Host " Log file: $LogFile" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

$Results = @()
$StartTime = Get-Date

foreach ($Sys in $Systems) {
    $SysStart = Get-Date
    Write-Host ""
    Write-Host "====================================================================" -ForegroundColor Cyan
    Write-Host " >>> [STARTING FORMAL 100 NS SIMULATION & ANALYSIS] System: $Sys" -ForegroundColor Cyan
    Write-Host "     Timestamp: $($SysStart.ToString('yyyy-MM-dd HH:mm:ss'))" -ForegroundColor Cyan
    Write-Host "====================================================================" -ForegroundColor Cyan

    $Status = "SUCCESS"
    "=== [FORMAL 100NS START] System: $Sys at $($SysStart.ToString('yyyy-MM-dd HH:mm:ss')) ===" | Out-File -FilePath $LogFile -Append -Encoding utf8

    & .\run_100ns_formal.ps1 -System $Sys 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $LogFile -Append
    if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
        $Status = "FAILED (Exit Code $LASTEXITCODE)"
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
    if ((Test-Path "..\md_ache\md.xtc") -and (Test-Path "..\md_$CompSys\md.xtc")) {
        Write-Host "  -> Comparing ache Monomer vs $CompSys Complex..." -ForegroundColor Yellow
        & .\plot_compare.ps1 -Protein "ache" -Complex $CompSys -OutDir "..\compare_ache_vs_$CompSys" 2>&1 | ForEach-Object { "$_" } | Tee-Object -FilePath $LogFile -Append
    }
}

$TotalTime = (Get-Date - $StartTime).ToString("dd\.hh\:mm\:ss")

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " 4-SYSTEM 100 NS FORMAL PRODUCTION MD & ANALYSIS MASTER REPORT" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
$Results | Format-Table -AutoSize
Write-Host "Total Pipeline Duration: $TotalTime" -ForegroundColor Green
Write-Host "Full execution log saved in: $LogFile" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
