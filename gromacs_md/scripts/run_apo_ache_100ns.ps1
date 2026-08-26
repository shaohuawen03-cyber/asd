# ============================================================
# TRUE APO CONTROL: re-run md_ache as STANDALONE AChE (no peptide)
#
# 背景 (v2.7.3 数据核查): 现有 md_ache 的拓扑里含有链 B = ALLLHRC 七肽,
# 即 md_ache 其实是复合物, 不是单体对照。本脚本用纯 AChE 单体
# (input/ache.pdb, 单链) 按与其余三个体系完全相同的 v1.0 协议
# (mdp/100ns, 100 ns) 重跑, 替换 md_ache。
#
# 安全设计:
#   * 现有 md_ache 会先整体改名备份为 md_ache_complex_backup
#     (只改名不复制, 瞬时完成, 不丢任何旧结果);
#   * 重跑前校验 input/ache.pdb 确实是单链 apo (无链 B), 否则中止;
#   * MD 结束后校验拓扑没有 topol_Protein_chain_B.itp, 否则报警;
#   * 分析用 v2.7.3 拓扑驱动的 index (单体模式, 无 Peptide 组),
#     contacts/bridging/氢键对比自动按 apo 处理。
#
# 用法 (在 gromacs_md\scripts 下):
#   .\run_apo_ache_100ns.ps1            # 正式 100 ns (mdp/100ns, 与其它三体系同协议)
#   .\run_apo_ache_100ns.ps1 -Testing   # 5000 步快速贯通验证 (mdp/test)
#   .\run_apo_ache_100ns.ps1 -OnlyAnalyze  # 只对已存在的 apo md_ache 做分析+出图
# ============================================================
param(
    [switch]$Testing,
    [switch]$OnlyAnalyze
)

try { chcp 65001 | Out-Null } catch {}
[Console]::OutputEncoding = [System.Text.Encoding]::UTF8
$OutputEncoding = [System.Text.Encoding]::UTF8
$ErrorActionPreference = "Continue"

$Here = $PSScriptRoot                                  # gromacs_md\scripts
$Root = Split-Path $Here -Parent                       # gromacs_md
$InputPdb = Join-Path $Root "input\ache.pdb"
$SrcPdb   = Join-Path $Root "input\alllhrc_complex.pdb"
$Work     = Join-Path $Root "md_ache"
$Backup   = Join-Path $Root "md_ache_complex_backup"

Write-Host "====================================================================" -ForegroundColor Cyan
Write-Host " TRUE APO CONTROL: re-run md_ache as standalone AChE (100 ns, v1.0 protocol)" -ForegroundColor Cyan
Write-Host "====================================================================" -ForegroundColor Cyan

$py = Get-Command python -ErrorAction SilentlyContinue
if (-not $py) { $py = Get-Command py -ErrorAction SilentlyContinue }
if (-not $py) { $py = Get-Command python3 -ErrorAction SilentlyContinue }

if (-not $OnlyAnalyze) {

    # ---------- 1. 准备并校验 apo 输入 PDB ----------
    if (Test-Path (Join-Path $Root "input\ache_complex.pdb")) {
        Write-Host "!! 警告: 存在 input\ache_complex.pdb —— run_all.sh 会优先使用它而不是 ache.pdb!" -ForegroundColor Red
        Write-Host "   若它不是纯 apo 单体, 请先移走/改名再运行本脚本。" -ForegroundColor Red
        exit 1
    }
    if (-not (Test-Path $InputPdb)) {
        if (-not (Test-Path $SrcPdb)) {
            Write-Host "ERROR: input\ache.pdb 与 input\alllhrc_complex.pdb 都不存在" -ForegroundColor Red
            exit 1
        }
        Write-Host ">> 从 $SrcPdb 提取纯 AChE 单体 (链 A) -> input\ache.pdb ..." -ForegroundColor Yellow
        Push-Location $Here
        try {
            & $py.Source ".\extract_ache_monomer.py" $SrcPdb $InputPdb
        } finally { Pop-Location }
        if ($LASTEXITCODE -ne 0 -or -not (Test-Path $InputPdb)) {
            Write-Host "ERROR: 提取 AChE 单体失败" -ForegroundColor Red
            exit 1
        }
    }

    Write-Host ">> 校验 input\ache.pdb 是否为单链 apo (无肽链) ..." -ForegroundColor Yellow
    $chains = @{}
    $resids = @{}
    foreach ($line in Get-Content $InputPdb -ErrorAction SilentlyContinue) {
        if ($line.Length -ge 26 -and ($line.StartsWith("ATOM") -or $line.StartsWith("HETATM"))) {
            $resname = $line.Substring(17, 3).Trim()
            if ($resname -in @("SOL", "HOH", "WAT", "NA", "CL")) { continue }
            $chain = $line.Substring(21, 1)
            $resseq = $line.Substring(22, 4).Trim()
            $chains[$chain] = $true
            $resids[$resseq] = $true
        }
    }
    $nChains = $chains.Count
    $nRes = $resids.Count
    Write-Host "   chains=$nChains  residues=$nRes" -ForegroundColor White
    if ($nChains -gt 1) {
        Write-Host "ERROR: input\ache.pdb 包含 $nChains 条链 ($($chains.Keys -join ', ')) —— 这不是单体对照!" -ForegroundColor Red
        Write-Host "       请先修正 input\ache.pdb (应只有链 A, 无链 B), 再运行本脚本。" -ForegroundColor Red
        exit 1
    }
    if ($nRes -lt 400 -or $nRes -gt 650) {
        Write-Host "ERROR: 残基数 $nRes 异常 (预期 AChE 单体约 530) —— 请检查 input\ache.pdb" -ForegroundColor Red
        exit 1
    }

    # ---------- 2. 备份现有 md_ache (若存在旧结果) ----------
    if (Test-Path $Work) {
        $hasOld = (Test-Path (Join-Path $Work "md.xtc")) -or (Test-Path (Join-Path $Work "md.tpr")) -or (Test-Path (Join-Path $Work "md_0_1.xtc"))
        if ($hasOld) {
            $dst = $Backup
            if (Test-Path $dst) {
                $dst = "$Backup-$(Get-Date -Format 'yyyyMMdd-HHmmss')"
            }
            Write-Host ">> 备份现有 md_ache (复合物结果) -> $dst" -ForegroundColor Yellow
            Move-Item -Path $Work -Destination $dst
            if (-not (Test-Path $dst)) {
                Write-Host "ERROR: 备份失败, 中止 (不会动原数据)" -ForegroundColor Red
                exit 1
            }
            Write-Host "   [OK] 旧复合物轨迹/结果已保留在 $dst (分析脚本不会读取它)" -ForegroundColor Green
        }
    }

    # ---------- 3. 跑 MD (与其它三体系相同的 v1.0 100ns 协议) ----------
    Write-Host ""
    Write-Host ">> 启动 MD: run_all.ps1 -System ache ($(if ($Testing) {'TESTING 5000步'} else {'正式 100 ns'})) ..." -ForegroundColor Cyan
    Write-Host "   (mdrun 期间请勿关闭窗口; 完成后自动进入分析)" -ForegroundColor Yellow
    Push-Location $Here
    try {
        if ($Testing) {
            & .\run_all.ps1 -System ache -Testing 2>&1 | ForEach-Object { "$_" }
        } else {
            & .\run_all.ps1 -System ache 2>&1 | ForEach-Object { "$_" }
        }
    } finally { Pop-Location }
    if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
        Write-Host "ERROR: run_all.ps1 退出码 $LASTEXITCODE" -ForegroundColor Red
        exit $LASTEXITCODE
    }

    # ---------- 4. MD 后校验: 拓扑必须无链 B ----------
    if (Test-Path (Join-Path $Work "topol_Protein_chain_B.itp")) {
        Write-Host "!! 警告: md_ache 拓扑里仍存在 topol_Protein_chain_B.itp —— 输入的 PDB 可能仍有第二条链!" -ForegroundColor Red
        Write-Host "   请检查 input\ache.pdb, 删除链 B 后重跑。" -ForegroundColor Red
        exit 1
    }
    Write-Host ">> [OK] 拓扑无链 B -> md_ache 是真正的 apo 单体对照" -ForegroundColor Green
}

# ---------- 5. 分析 + 出图 (v2.7.3 拓扑驱动 index) ----------
Write-Host ""
Write-Host ">> 分析 + 出图: run_analysis.ps1 -System ache (monomer 模式) ..." -ForegroundColor Cyan
Push-Location $Here
try {
    & .\run_analysis.ps1 -System ache 2>&1 | ForEach-Object { "$_" }
} finally { Pop-Location }
if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
    Write-Host "ERROR: run_analysis.ps1 退出码 $LASTEXITCODE" -ForegroundColor Red
    exit $LASTEXITCODE
}

# ---------- 6. 刷新统一图 + 对照图 (ache 自动按 apo 处理) ----------
Write-Host ""
Write-Host ">> 刷新四体系统一图与对照图 (run_unified_replot.ps1) ..." -ForegroundColor Cyan
Push-Location $Here
try {
    & .\run_unified_replot.ps1 2>&1 | ForEach-Object { "$_" }
} finally { Pop-Location }

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " DONE." -ForegroundColor Green
Write-Host "   md_ache                      = 新 apo 单体对照 (100 ns)" -ForegroundColor Green
Write-Host "   md_ache_complex_backup*      = 旧复合物结果 (未删除)" -ForegroundColor Yellow
Write-Host " 检查:" -ForegroundColor Green
Write-Host "   md_ache\figures\fig0_summary_all.png  -> 'AChE MD summary (apo control...)'" -ForegroundColor Green
Write-Host "   compare_ache_vs_*\fig_compare.png     -> F 面板只有复合物的 AChE-Peptide 曲线, ache 不出现" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
