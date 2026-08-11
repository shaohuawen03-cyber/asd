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
Write-Host " 启动四套体系 (alllhrc, fllhttr, ylsllqr, ache) 5000 步全流程贯通验证" -ForegroundColor Green
Write-Host " 日志保存路径: $LogFile" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

$Results = @()

foreach ($Sys in $Systems) {
    Write-Host ""
    Write-Host "====================================================================" -ForegroundColor Cyan
    Write-Host " >>> 正在运行全流程贯通验证体系: $Sys (从最底层建系到5000步MD到绘图)" -ForegroundColor Cyan
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
Write-Host " 4-SYSTEM FULL WORKFLOW TEST REPORT (四体系贯通验证结果汇总表)" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
$Results | Format-Table -AutoSize
Write-Host "完整运行日志已保存至: $LogFile" -ForegroundColor Green

$Errors = Select-String -Path $LogFile -Pattern "Fatal error|Segmentation fault|command not found|No such file|ModuleNotFoundError|Syntax error" -ErrorAction SilentlyContinue
if ($Errors) {
    Write-Host ""
    Write-Host ">>> [在日志 test_four_systems.log 中检测到报错信息] :" -ForegroundColor Red
    $Errors | ForEach-Object { Write-Host "  $_.Line" -ForegroundColor Yellow }
    Write-Host "`n请复制上面黄色的具体报错行发给我！" -ForegroundColor Red
} else {
    Write-Host ""
    Write-Host ">>> [大获成功！4 套体系 100% 顺畅通关 - 日志中零 ERROR！]" -ForegroundColor Green
}
Write-Host "====================================================================" -ForegroundColor Green
