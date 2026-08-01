# ============================================================
# AChE (4ey6) - beta-淀粉样肽(Aβ)复合物 分子动力学模拟 分析全流程 (PowerShell 原生版)
# 复现论文《乙酰胆碱酯酶-β-淀粉样肽复合物的分子动力学模拟》第 3 部分
#
# 用法:  .\run_analysis.ps1 -System alllhrc -Testing
#        .\run_analysis.ps1 -System alllhrc
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System,

    [switch]$Testing
)

$ErrorActionPreference = "Stop"

# 自动检测 gmx
$GMX = "gmx"
if (Get-Command "gmx.exe" -ErrorAction SilentlyContinue) {
    $GMX = "gmx.exe"
}

# 自动检测 python (在 Windows Anaconda 中执行程序名为主解释器 python)
$PY = "python"
if (Get-Command "python.exe" -ErrorAction SilentlyContinue) {
    $PY = "python.exe"
}

$WorkDir = "..\md_$System"
if (-not (Test-Path "$WorkDir\md.xtc")) {
    Write-Host "!!! 错误: 在 $WorkDir 目录下未找到产物轨迹文件 md.xtc！" -ForegroundColor Red
    Write-Host "!!! 请先执行 MD 模拟生成该轨迹。" -ForegroundColor Red
    exit 1
}

if ($Testing -or $env:TESTING -eq "1") {
    Write-Host ">> [分析测试模式] 适应 100 步测试轨迹的自适应统计参数" -ForegroundColor Yellow
    $env:TESTING = "1"
} else {
    Write-Host ">> [分析正式模式] 针对 1000 ns 生产轨迹执行完整计算" -ForegroundColor Cyan
    $env:TESTING = "0"
}

Push-Location $WorkDir
try {
    # 取 Scripts 分析脚本的绝对路径 (标准 Windows 路径格式)
    $ScriptsDir = (Resolve-Path "..\scripts\analysis").Path

    Write-Host "========== 开始对体系 $System 执行论文结果分析 (PowerShell) ========" -ForegroundColor Green
    Write-Host ">> 采用 Python 解释器: $PY" -ForegroundColor Green

    Write-Host "[0/8] 建立分析分组索引 (0_make_index.sh) ..."
    & bash "$ScriptsDir/0_make_index.sh"

    Write-Host "[1/8] 计算骨架 RMSD / RMSF (1_rmsd_rmsf.sh) ..."
    & bash "$ScriptsDir/1_rmsd_rmsf.sh"

    Write-Host "[2/8] 计算 Aβ 围绕 AChE 的径向分布函数 RDF (2_rdf.sh) ..."
    & bash "$ScriptsDir/2_rdf.sh"

    Write-Host "[3/8] 计算复合物 SASA (3_sasa.sh) ..."
    & bash "$ScriptsDir/3_sasa.sh"

    Write-Host "[4/8] 计算 Aβ 肽二级结构演变 (4_secondary_structure.sh) ..."
    & bash "$ScriptsDir/4_secondary_structure.sh"

    Write-Host "[5/8] 统计间/内氢键分布 (5_hbond.sh) ..."
    & bash "$ScriptsDir/5_hbond.sh"

    Write-Host "[6/8] 统计天然与非天然相互作用接触 (contacts.py) ..."
    & $PY "$ScriptsDir\contacts.py" -t md.tpr -f md.xtc

    Write-Host "[7/8] 统计水介导桥连相互作用 (bridging_waters.py) ..."
    & $PY "$ScriptsDir\bridging_waters.py" -t md.tpr -f md.xtc

    Write-Host "[8/8] 批量绘制论文出版级图表 (SVG / PNG / PDF) ..."
    & $PY "$ScriptsDir\plot_all.py" --dir . --out ./figures

    Write-Host "========== 体系 $System 分析流程全部完成！ ==========" -ForegroundColor Green
    Write-Host "生成分析图表一览 (保存在 ./figures/ 下):"
    Write-Host "  - fig1_rmsd_rmsf.{svg,png,pdf}            => 图 1 (RMSD / RMSF)"
    Write-Host "  - fig2_rdf.{svg,png,pdf}                  => 图 2 (径向分布函数 RDF)"
    Write-Host "  - fig3_sasa.{svg,png,pdf}                 => 图 3 (溶剂可及表面积 SASA)"
    Write-Host "  - fig4_secondary_structure.{svg,png,pdf}  => 图 4 (二级结构倾向)"
    Write-Host "  - fig5_contacts.{svg,png,pdf}             => 图 5 (非天然残基对接触)"
    Write-Host "  - fig6_bridging_waters.{svg,png,pdf}      => 图 6 (水介导桥连)"
    Write-Host "  - fig_hbonds.{svg,png,pdf}                => 论文 3.3 节 (氢键数量曲线)"
    Write-Host "  - fig0_summary_all.{svg,png,pdf}          => 综合 2x3 六格汇总对比主图"
    Write-Host "  - summary_metrics.csv & wide.csv          => 统计指标汇总表"
} finally {
    Pop-Location
}
