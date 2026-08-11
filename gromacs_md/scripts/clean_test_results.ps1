# ============================================================
# Clean Up All Previous Test Results and Temporary Files (PowerShell)
#
# Usage:  .\clean_test_results.ps1
# ============================================================
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::InputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

$DirsToClean = @("..\md_alllhrc", "..\md_fllhttr", "..\md_ylsllqr", "..\md_ache", "..\compare_figures")
$FilesToClean = @(".\test_four_systems.log")

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " Cleaning up all previous test directories and log files..." -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

foreach ($Dir in $DirsToClean) {
    if (Test-Path $Dir) {
        Remove-Item -Path $Dir -Recurse -Force -ErrorAction SilentlyContinue
        Write-Host "  [DELETED DIRECTORY] $Dir" -ForegroundColor Yellow
    }
}

foreach ($File in $FilesToClean) {
    if (Test-Path $File) {
        Remove-Item -Path $File -Force -ErrorAction SilentlyContinue
        Write-Host "  [DELETED FILE] $File" -ForegroundColor Yellow
    }
}

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " SUCCESS! Workspace is completely clean and ready for formal 100 ns simulations." -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
