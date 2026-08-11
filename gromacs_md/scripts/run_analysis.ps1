# ============================================================
# AChE - A-beta Complex Trajectory Analysis & Figure Generation (PowerShell)
#
# Usage:  .\run_analysis.ps1 -System alllhrc
#         .\run_analysis.ps1 -System alllhrc -Testing
#         .\run_analysis.ps1 -System alllhrc -OnlyPlot
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System,

    [switch]$Testing,
    [switch]$OnlyPlot
)

[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
[Console]::InputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

$GMX = "gmx"
if (Get-Command "gmx.exe" -ErrorAction SilentlyContinue) {
    $GMX = "gmx.exe"
}

$PY = "python"
if (Get-Command "python.exe" -ErrorAction SilentlyContinue) {
    $PY = "python.exe"
}

$WorkDir = "..\md_$System"
if (-not (Test-Path "$WorkDir\md.xtc")) {
    Write-Host "!!! 错误: 在 $WorkDir 目录下未找到产物轨迹文件 md.xtc！" -ForegroundColor Red
    Write-Host "!!! 请先运行动力学模拟。" -ForegroundColor Red
    exit 1
}

if ($Testing -or $env:TESTING -eq "1") {
    Write-Host ">> [分析测试模式] 针对 5000 步测试轨迹自动适配计算" -ForegroundColor Yellow
    $env:TESTING = "1"
} else {
    Write-Host ">> [分析正式模式] 针对 100 ns 生产轨迹执行完整统计" -ForegroundColor Cyan
    $env:TESTING = "0"
}

Push-Location $WorkDir
try {
    Write-Host "========== 开始对体系 $System 执行论文分析计算 ========" -ForegroundColor Green
    Write-Host ">> 采用 Python 解释器: $PY" -ForegroundColor Green

    if (-not $OnlyPlot) {
        Write-Host "[0/8] 建立分析分组索引 (0_make_index.sh) ..."
        & bash "../scripts/analysis/0_make_index.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[1/8] 计算骨架 RMSD / RMSF (1_rmsd_rmsf.sh) ..."
        & bash "../scripts/analysis/1_rmsd_rmsf.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[2/8] 计算肽围绕 AChE 质心径向分布函数 RDF (2_rdf.sh) ..."
        & bash "../scripts/analysis/2_rdf.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[3/8] 计算溶剂可及表面积 SASA (3_sasa.sh) ..."
        & bash "../scripts/analysis/3_sasa.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[4/8] 计算肽二级结构分布演变 (4_secondary_structure.sh) ..."
        & bash "../scripts/analysis/4_secondary_structure.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[5/8] 统计界面与内部氢键分布 (5_hbond.sh) ..."
        & bash "../scripts/analysis/5_hbond.sh" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[6/8] 统计天然与非天然相互作用接触 (contacts.py) ..."
        & $PY "..\scripts\analysis\contacts.py" -t md.tpr -f md.xtc 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        Write-Host "[7/8] 统计水介导桥连相互作用 (bridging_waters.py) ..."
        & $PY "..\scripts\analysis\bridging_waters.py" -t md.tpr -f md.xtc 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    } else {
        Write-Host ">> [OnlyPlot 模式] 跳过 0-7 步 GROMACS 计算，直接刷新全部图像 ..." -ForegroundColor Yellow
    }

    Write-Host "[8/8] 批量绘制论文出版级图表 (SVG / PNG / PDF) ..."
    & $PY "..\scripts\analysis\plot_all.py" --dir . --out ./figures 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    Write-Host "========== 体系 $System 统计与图表生成全部完成！ ==========" -ForegroundColor Green
    Write-Host "生成分析图表一览 (保存在 ./figures/ 下):"
    Write-Host "  - fig1_rmsd_rmsf.{svg,png,pdf}            => 图 1 (RMSD / RMSF)"
    Write-Host "  - fig2_rdf.{svg,png,pdf}                  => 图 2 (径向分布函数 RDF)"
    Write-Host "  - fig3_sasa.{svg,png,pdf}                 => 图 3 (溶剂可及表面积 SASA)"
    Write-Host "  - fig4_secondary_structure.{svg,png,pdf}  => 图 4 (二级结构倾向)"
    Write-Host "  - fig5_contacts.{svg,png,pdf}             => 图 5 (非天然残基对接触)"
    Write-Host "  - fig6_bridging_waters.{svg,png,pdf}      => 图 6 (水介导桥连)"
    Write-Host "  - fig_hbonds.{svg,png,pdf}                => 论文 3.3 节 (氢键数量曲线)"
    Write-Host "  - fig0_summary_all.{svg,png,pdf}          => 综合 2x3 六格组合大总图"
    Write-Host "  - summary_metrics.csv & wide.csv          => 统计指标汇总表"
} finally {
    Pop-Location
}
