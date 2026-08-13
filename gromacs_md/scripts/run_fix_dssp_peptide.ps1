# ============================================================
# 只重跑「有问题的部分」: 肽 DSSP + 肽 RMSD 三段相位 + 重新出图
# 不重跑 RDF / SASA / 氢键 / 接触 / 桥连水, 也不动正在跑的 mdrun
#
# Usage:  .\run_fix_dssp_peptide.ps1 -System alllhrc
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System
)

[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::InputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

$PY = "python"
if (Get-Command "python.exe" -ErrorAction SilentlyContinue) {
    $PY = "python.exe"
}

$WorkDir = "..\md_$System"
if (-not (Test-Path $WorkDir)) {
    Write-Host "!!! 错误: 找不到工作目录 $WorkDir" -ForegroundColor Red
    exit 1
}
if (-not (Test-Path "$WorkDir\md_0_1.xtc") -and -not (Test-Path "$WorkDir\md.xtc") -and -not (Test-Path "$WorkDir\md_fit.xtc")) {
    Write-Host "!!! 错误: $WorkDir 下没有轨迹 (md_fit.xtc / md_0_1.xtc / md.xtc)" -ForegroundColor Red
    exit 1
}

Write-Host "========== 只修复 DSSP + 肽三段 RMSD: $System ==========" -ForegroundColor Green
Write-Host "不会重跑 RDF/SASA/Hbond/Contacts, 也不会中断后台 mdrun" -ForegroundColor Yellow

Push-Location $WorkDir
try {
    if (-not (Test-Path "index.ndx")) {
        Write-Host "[0] index.ndx 缺失, 补建分析分组 (0_make_index.sh) ..."
        & bash "../scripts/analysis/0_make_index.sh" 2>&1 | ForEach-Object { "$_" }
    } else {
        Write-Host "[0] 已有 index.ndx, 跳过 make_ndx / trjconv" -ForegroundColor Yellow
    }

    Write-Host "[1/3] 肽 DSSP 逐帧重建 (4_secondary_structure.sh -> compute_peptide_ss.py) ..."
    & bash "../scripts/analysis/4_secondary_structure.sh" 2>&1 | ForEach-Object { "$_" }

    Write-Host "[2/3] 肽 RMSD 三段相位 + 可选 AChE-fitted 配体 RMSD ..."
    & $PY "..\scripts\analysis\analyze_peptide_phases.py" -d . 2>&1 | ForEach-Object { "$_" }

    Write-Host "[3/3] 仅重绘出版图 (plot_all.py, 读取已有 xvg, 不重算其它分析) ..."
    & $PY "..\scripts\analysis\plot_all.py" --dir . --out ./figures 2>&1 | ForEach-Object { "$_" }
} finally {
    Pop-Location
}

Write-Host "========== DSSP + 肽三段 修复完成: md_$System\figures\ ==========" -ForegroundColor Green
Write-Host "请检查:"
Write-Host "  fig1_rmsd_rmsf.*            (P1/P2/P3 色带标注三段)"
Write-Host "  fig4_secondary_structure.*  (逐帧 DSSP + 逐残基热图)"
Write-Host "  fig0_summary_all.*          (E 面板不再写 N/A for Monomer)"
Write-Host "  ss_pep_summary.txt          (DSSP 文字报告)"
Write-Host "  peptide_rmsd_phases.txt     (三段均值表)"
