# ============================================================
# AChE & A-beta Four Systems 5000-Step Automation Suite (PowerShell)
#
# Usage:  .\test_four_systems.ps1
# ============================================================
$Systems = @("alllhrc", "fllhttr", "ylsllqr", "ache")

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " Starting 4-System (alllhrc, fllhttr, ylsllqr, ache) 5000-Step MD Suite" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

foreach ($Sys in $Systems) {
    Write-Host ""
    Write-Host "====================================================================" -ForegroundColor Cyan
    Write-Host " Processing System: $Sys (5000-Step Validation)" -ForegroundColor Cyan
    Write-Host "====================================================================" -ForegroundColor Cyan
    
    .\run_pipeline_all.ps1 -System $Sys -Testing
}

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " SUCCESS! All 4 systems completed MD simulation and analysis plots!" -ForegroundColor Green
Write-Host " Check generated figures in: md_<System>/figures/" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
