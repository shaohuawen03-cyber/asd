# ============================================================
# AChE (4ey6) - beta-淀粉样肽(Aβ)复合物 分子动力学模拟 全流程 (PowerShell 原生脚本)
#
# 用法:  .\run_all.ps1 -System alllhrc -Testing
#        .\run_all.ps1 -System alllhrc
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System,

    [switch]$Testing
)

$ErrorActionPreference = "Stop"

if ($Testing -or $env:TESTING -eq "1") {
    Write-Host ">> [测试模式] 运行 100 步测试流程" -ForegroundColor Yellow
    $env:TESTING = "1"
} else {
    Write-Host ">> [正式模式] 运行 1000 ns 产物流程" -ForegroundColor Cyan
    $env:TESTING = "0"
}

# 自动调用 Bash 执行核心流程
& bash ".\run_all.sh" $System
