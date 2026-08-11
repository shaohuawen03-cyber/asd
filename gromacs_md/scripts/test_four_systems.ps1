# ============================================================
# AChE & Aβ肽 四套体系 (三个对接复合物 + AChE 单体) 10,000 步测试自动化套件
#
# 用法:  .\test_four_systems.ps1
# ============================================================
$Systems = @("alllhrc", "fllhttr", "ylsllqr", "ache")

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " 启动四体系 (alllhrc, fllhttr, ylsllqr, ache) 10000 步全自动化动力学与图表测试" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

foreach ($Sys in $Systems) {
    Write-Host "`n====================================================================" -ForegroundColor Cyan
    Write-Host " >>> 正在处理体系: $Sys (10000 步验证)" -ForegroundColor Cyan
    Write-Host "====================================================================" -ForegroundColor Cyan
    
    .\run_pipeline_all.ps1 -System $Sys -Testing
}

Write-Host "`n====================================================================" -ForegroundColor Green
Write-Host " 恭喜！四套体系 (三个复合物 + 单体 AChE) 10000 步 MD 模拟及论文分析图表全数完成！" -ForegroundColor Green
Write-Host "  - 查看各体系出图: md_<体系名>\figures\" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
