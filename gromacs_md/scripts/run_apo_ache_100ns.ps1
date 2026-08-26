# ============================================================
# TRUE APO CONTROL: re-run md_ache as STANDALONE AChE (no peptide)
#
# Background (v2.7.3 data check): the current md_ache topology contains
# chain B = the ALLLHRC 7-mer peptide, i.e. md_ache is actually a complex,
# not a protein-only control. This script re-runs md_ache as a pure AChE
# monomer (input/ache.pdb, single chain) with the SAME v1.0 protocol as the
# other three systems (mdp/100ns, 100 ns) and replaces md_ache.
#
# Safety design:
#   * the existing md_ache is first RENAMED to md_ache_complex_backup
#     (rename only, instant, nothing is deleted);
#   * input/ache.pdb is verified to be a single-chain apo PDB (no chain B)
#     BEFORE anything runs; otherwise the script aborts;
#   * after MD the topology is checked for topol_Protein_chain_B.itp
#     (must NOT exist);
#   * analysis uses the v2.7.3 topology-driven monomer index (no Peptide
#     group), so contacts/bridging/H-bond comparisons treat ache as apo
#     automatically.
#
# NOTE: this file is intentionally ASCII-only (English). Windows PowerShell
# 5.1 reads BOM-less .ps1 files with the system ANSI codepage (GBK on
# Chinese Windows), which garbles UTF-8 Chinese text and breaks parsing.
#
# Usage (from gromacs_md\scripts):
#   .\run_apo_ache_100ns.ps1            # formal 100 ns (mdp/100ns)
#   .\run_apo_ache_100ns.ps1 -Testing   # 5000-step quick validation (mdp/test)
#   .\run_apo_ache_100ns.ps1 -OnlyAnalyze  # analyze+plot existing apo md_ache
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

function Get-PdbChainInfo {
    param([string]$Path)
    $chains = @{}
    $resids = @{}
    foreach ($line in Get-Content $Path -ErrorAction SilentlyContinue) {
        if ($line.Length -ge 26 -and ($line.StartsWith("ATOM") -or $line.StartsWith("HETATM"))) {
            $resname = $line.Substring(17, 3).Trim()
            if ($resname -in @("SOL", "HOH", "WAT", "NA", "CL")) { continue }
            $chain = $line.Substring(21, 1).Trim()
            $resseq = $line.Substring(22, 4).Trim()
            $chains[$chain] = $true
            $resids[$resseq] = $true
        }
    }
    return @{ Chains = $chains; Residues = $resids }
}

function Test-ApoPdb {
    param([string]$Path)
    $info = Get-PdbChainInfo $Path
    $ok = ($info.Chains.Count -eq 1) -and ($info.Residues.Count -ge 400) -and ($info.Residues.Count -le 650)
    return @{ Ok = $ok; Chains = $info.Chains.Count; Residues = $info.Residues.Count }
}

if (-not $OnlyAnalyze) {

    # ---------- 1. prepare and validate the apo input PDB ----------
    if (Test-Path (Join-Path $Root "input\ache_complex.pdb")) {
        Write-Host "!! WARNING: input\ache_complex.pdb exists - run_all.sh would use it instead of ache.pdb!" -ForegroundColor Red
        Write-Host "   If it is not a pure apo monomer, move/rename it before running this script." -ForegroundColor Red
        exit 1
    }

    $apoOk = $false
    if (Test-Path $InputPdb) {
        $chk = Test-ApoPdb $InputPdb
        if ($chk.Ok) {
            $apoOk = $true
            Write-Host ">> [OK] input\ache.pdb: single chain, $($chk.Residues) residues" -ForegroundColor Green
        } else {
            Write-Host "!! input\ache.pdb is NOT a single-chain apo PDB (chains=$($chk.Chains), residues=$($chk.Residues))." -ForegroundColor Yellow
            Write-Host "   It still contains a second chain (probably the peptide). Rebuilding it from chain A only." -ForegroundColor Yellow
        }
    } else {
        Write-Host ">> input\ache.pdb not found - will extract chain A from alllhrc_complex.pdb." -ForegroundColor Yellow
    }

    if (-not $apoOk) {
        if (-not (Test-Path $SrcPdb)) {
            Write-Host "ERROR: input\alllhrc_complex.pdb not found - cannot rebuild the apo PDB." -ForegroundColor Red
            Write-Host "       Provide a single-chain AChE PDB at input\ache.pdb (no chain B), or restore the source complex PDB, then re-run." -ForegroundColor Red
            exit 1
        }
        if (-not $py) {
            Write-Host "ERROR: python not found on PATH - needed to extract chain A." -ForegroundColor Red
            exit 1
        }
        if (Test-Path $InputPdb) {
            $bakPdb = "${InputPdb}.bak-$(Get-Date -Format 'yyyyMMdd-HHmmss')"
            Move-Item -Path $InputPdb -Destination $bakPdb
            Write-Host "   [BACKUP] old input\ache.pdb -> $bakPdb" -ForegroundColor Yellow
        }
        Write-Host ">> extracting chain A only from $SrcPdb -> input\ache.pdb ..." -ForegroundColor Cyan
        Push-Location $Here
        try {
            & $py.Source ".\extract_ache_monomer.py" $SrcPdb $InputPdb
        } finally { Pop-Location }
        if ($LASTEXITCODE -ne 0 -or -not (Test-Path $InputPdb)) {
            Write-Host "ERROR: extraction of the AChE monomer failed." -ForegroundColor Red
            exit 1
        }
        $chk = Test-ApoPdb $InputPdb
        if (-not $chk.Ok) {
            Write-Host "ERROR: rebuilt input\ache.pdb is still not a single-chain apo PDB (chains=$($chk.Chains), residues=$($chk.Residues))." -ForegroundColor Red
            Write-Host "       Check ${SrcPdb}: chain A should be AChE, chain B the peptide." -ForegroundColor Red
            exit 1
        }
        Write-Host ">> [OK] rebuilt input\ache.pdb: single chain, $($chk.Residues) residues" -ForegroundColor Green
    }

    # ---------- 2. back up the existing md_ache (complex results) ----------
    if (Test-Path $Work) {
        $hasOld = (Test-Path (Join-Path $Work "md.xtc")) -or (Test-Path (Join-Path $Work "md.tpr")) -or (Test-Path (Join-Path $Work "md_0_1.xtc"))
        if ($hasOld) {
            $dst = $Backup
            if (Test-Path $dst) {
                $dst = "$Backup-$(Get-Date -Format 'yyyyMMdd-HHmmss')"
            }
            Write-Host ">> backing up current md_ache (complex results) -> $dst" -ForegroundColor Yellow
            Move-Item -Path $Work -Destination $dst
            if (-not (Test-Path $dst)) {
                Write-Host "ERROR: backup failed - aborting (original data untouched)" -ForegroundColor Red
                exit 1
            }
            Write-Host "   [OK] old complex trajectory/results kept in $dst (analysis scripts ignore it)" -ForegroundColor Green
        }
    }

    # ---------- 3. run MD (same v1.0 100 ns protocol as the other systems) ----------
    Write-Host ""
    if ($Testing) {
        Write-Host ">> starting MD: run_all.ps1 -System ache -Testing (5000-step validation) ..." -ForegroundColor Cyan
    } else {
        Write-Host ">> starting MD: run_all.ps1 -System ache (formal 100 ns) ..." -ForegroundColor Cyan
    }
    Write-Host "   (keep this window open while mdrun runs; analysis starts automatically afterwards)" -ForegroundColor Yellow
    Push-Location $Here
    try {
        if ($Testing) {
            & .\run_all.ps1 -System ache -Testing 2>&1 | ForEach-Object { "$_" }
        } else {
            & .\run_all.ps1 -System ache 2>&1 | ForEach-Object { "$_" }
        }
    } finally { Pop-Location }
    if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
        Write-Host "ERROR: run_all.ps1 exited with code $LASTEXITCODE" -ForegroundColor Red
        exit $LASTEXITCODE
    }

    # ---------- 4. post-MD check: topology must have NO chain B ----------
    if (Test-Path (Join-Path $Work "topol_Protein_chain_B.itp")) {
        Write-Host "!! WARNING: md_ache topology still has topol_Protein_chain_B.itp - the input PDB probably still has a second chain!" -ForegroundColor Red
        Write-Host "   Fix input\ache.pdb (remove chain B), then re-run." -ForegroundColor Red
        exit 1
    }
    Write-Host ">> [OK] no chain B in topology -> md_ache is a true apo monomer control" -ForegroundColor Green
}

# ---------- 5. analysis + figures (v2.7.3 topology-driven index) ----------
Write-Host ""
Write-Host ">> analysis + plotting: run_analysis.ps1 -System ache (monomer mode) ..." -ForegroundColor Cyan
Push-Location $Here
try {
    & .\run_analysis.ps1 -System ache 2>&1 | ForEach-Object { "$_" }
} finally { Pop-Location }
if ($LASTEXITCODE -ne 0 -and $null -ne $LASTEXITCODE) {
    Write-Host "ERROR: run_analysis.ps1 exited with code $LASTEXITCODE" -ForegroundColor Red
    exit $LASTEXITCODE
}

# ---------- 6. refresh unified figures + compare folders (ache auto-handled as apo) ----------
Write-Host ""
Write-Host ">> refreshing unified figures + compare folders (run_unified_replot.ps1) ..." -ForegroundColor Cyan
Push-Location $Here
try {
    & .\run_unified_replot.ps1 2>&1 | ForEach-Object { "$_" }
} finally { Pop-Location }

Write-Host ""
Write-Host "====================================================================" -ForegroundColor Green
Write-Host " DONE." -ForegroundColor Green
Write-Host "   md_ache                   = new apo monomer control (100 ns)" -ForegroundColor Green
Write-Host "   md_ache_complex_backup*   = old complex results (kept, not deleted)" -ForegroundColor Yellow
Write-Host " Check:" -ForegroundColor Green
Write-Host "   md_ache\figures\fig0_summary_all.png  -> 'AChE MD summary (apo control...)'" -ForegroundColor Green
Write-Host "   compare_ache_vs_*\fig_compare.png     -> panel F has only the complex's AChE-Peptide curve (no ache)" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
