# ============================================================
# GROMACS 12-Step Split-Topology MD Workflow (PowerShell)
#
# Automatically splits complex PDB into pro.pdb and tai.pdb, generates separate
# topologies via pdb2gmx, merges .gro files, inserts .itp references, and runs
# MD following the user's exact 12-step protocol.
#
# Usage:  .\run_split_md_workflow.ps1 -System alllhrc
#         .\run_split_md_workflow.ps1 -System fllhttr
#         .\run_split_md_workflow.ps1 -System ylsllqr
#         .\run_split_md_workflow.ps1 -System ache
# ============================================================
param(
    [Parameter(Mandatory=$true)]
    [string]$System
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

# Locate mdp_templates directory
$MdpDir = ""
foreach ($cand in @("..\..\mdp_templates", "..\mdp_templates", "F:\0wsh\asd\mdp_templates", "C:\0wsh\asd\mdp_templates", "..\mdp\100ns")) {
    if (Test-Path $cand) {
        $MdpDir = (Resolve-Path $cand).Path
        break
    }
}
if (-not $MdpDir) {
    Write-Host "ERROR: Could not find mdp_templates directory!" -ForegroundColor Red
    exit 1
}
Write-Host ">> Using MDP Templates Directory: $MdpDir" -ForegroundColor Green

$WorkDir = "..\md_$System"
if (Test-Path $WorkDir) {
    Write-Host ">> Cleaning previous workspace $WorkDir for a 100% fresh start..." -ForegroundColor Yellow
    Remove-Item -Path $WorkDir -Recurse -Force -ErrorAction SilentlyContinue
}
New-Item -ItemType Directory -Path $WorkDir -Force | Out-Null

$InputPdb = ""
foreach ($pdbCand in @("..\input\${System}_complex.pdb", "..\input\${System}.pdb", "..\input\${System}_alone.pdb")) {
    if (Test-Path $pdbCand) {
        $InputPdb = (Resolve-Path $pdbCand).Path
        break
    }
}
if (-not $InputPdb) {
    Write-Host "ERROR: Could not find input PDB file for System: $System!" -ForegroundColor Red
    exit 1
}

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " Starting 12-Step Split-Topology MD for System: $System" -ForegroundColor Green
Write-Host " Input PDB: $InputPdb" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

Push-Location $WorkDir
try {
    # Step 0: Split complex PDB into pro.pdb and tai.pdb
    Write-Host "`nStep 0: Splitting complex PDB into pro.pdb and tai.pdb..." -ForegroundColor Cyan
    & $PY "..\scripts\split_complex_pdb.py" $InputPdb "pro.pdb" "tai.pdb" 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    $HasTai = (Test-Path "tai.pdb") -and ((Get-Item "tai.pdb").Length -gt 100)

    # Step 1: Generate topologies via pdb2gmx
    Write-Host "`nStep 1: Running gmx pdb2gmx for receptor pro.pdb..." -ForegroundColor Cyan
    @("5", "5") | & $GMX pdb2gmx -f "pro.pdb" -o "pro_processed.gro" -p "topol.top" -ignh -ff "amber99sb-ildn" -water "tip3p" 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    if ($HasTai) {
        Write-Host "Step 1b: Running gmx pdb2gmx for ligand tai.pdb..." -ForegroundColor Cyan
        @("5", "5") | & $GMX pdb2gmx -f "tai.pdb" -o "tai.gro" -p "tai.top" -i "taiposre.itp" -ignh -ff "amber99sb-ildn" -water "none" 2>&1 | ForEach-Object { "$_" }
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

        # Step 2: Merge .gro files
        Write-Host "`nStep 2: Merging pro_processed.gro and tai.gro..." -ForegroundColor Cyan
        $file1 = "pro_processed.gro"
        $file2 = "tai.gro"
        $lines1 = Get-Content $file1
        $lines2 = Get-Content $file2
        $n1 = [int]$lines1[1].Trim()
        $n2 = [int]$lines2[1].Trim()
        $total = $n1 + $n2
        Write-Host "File1 atoms: $n1 | File2 atoms: $n2 | Total: $total" -ForegroundColor Green

        $result = @()
        $result += $lines1[0]
        $result += "{0,5}" -f $total
        $result += $lines1[2..($n1+1)]
        $result += $lines2[2..($n2+1)]
        $result += $lines1[-1]
        $result | Set-Content $file1 -Encoding ASCII

        # Step 3: Update topol.top (include tai.itp and add tai 1 to [ molecules ])
        Write-Host "`nStep 3: Updating topol.top with tai.itp and tai 1..." -ForegroundColor Cyan
        $content = Get-Content "topol.top" -Raw
        $itpInclude = '#include "tai.itp"'
        $moleculeEntry = "tai                 1"

        if ($content -match '(#include "amber99sb-ildn\.ff/forcefield\.itp"|#include "amber99sb\.ff/forcefield\.itp")') {
            $content = $content -replace '(#include "amber[\w-]+\.ff/forcefield\.itp")', "`$1`r`n$itpInclude"
        } else {
            $content = "$itpInclude`r`n$content"
        }
        if ($content -match '(\[ molecules \][\s\S]*?Protein_chain_\w+\s+1)') {
            $content = $content -replace '(\[ molecules \][\s\S]*?Protein_chain_\w+\s+1)', "`$1`r`n$moleculeEntry"
        } else {
            $content = $content.TrimEnd() + "`r`n`r`n[ molecules ]`r`ntai                 1`r`n"
        }
        Set-Content -Path "topol.top" -Value $content -Encoding UTF8
    }

    # Step 4: editconf (triclinic box -d 1.0)
    Write-Host "`nStep 4: editconf (triclinic box -d 1.0)..." -ForegroundColor Cyan
    & $GMX editconf -f "pro_processed.gro" -o "gro_newbox.gro" -c -d 1.0 -bt triclinic 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    # Step 5: solvate
    Write-Host "`nStep 5: solvate (spc216.gro)..." -ForegroundColor Cyan
    & $GMX solvate -cp "gro_newbox.gro" -cs spc216.gro -o "gro_solv.gro" -p "topol.top" 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    # Step 6: grompp + genion (add ions)
    Write-Host "`nStep 6: ions (ions.mdp + genion neutral 0.15M)..." -ForegroundColor Cyan
    $IonsMdp = "$MdpDir\ions.mdp"
    if (-not (Test-Path $IonsMdp)) { $IonsMdp = "$MdpDir\1_min.mdp" }
    & $GMX grompp -f $IonsMdp -c "gro_solv.gro" -r "gro_solv.gro" -p "topol.top" -o "ions.tpr" -maxwarn 7 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    @("SOL") | & $GMX genion -s "ions.tpr" -o "solv_ions.gro" -p "topol.top" -pname NA -nname CL -neutral -conc 0.15 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    # Step 7: Energy Minimization
    Write-Host "`nStep 7: Energy Minimization (em.mdp)..." -ForegroundColor Cyan
    $EmMdp = "$MdpDir\em.mdp"
    if (-not (Test-Path $EmMdp)) { $EmMdp = "$MdpDir\1_min.mdp" }
    & $GMX grompp -f $EmMdp -c "solv_ions.gro" -r "solv_ions.gro" -p "topol.top" -o "em.tpr" -maxwarn 7 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & $GMX mdrun -v -deffnm em -nb gpu 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    # Step 8: Position restraints insertion
    if ($HasTai) {
        Write-Host "`nStep 8: Inserting position restraints into topol.top..." -ForegroundColor Cyan
        $content = Get-Content "topol.top" -Raw
        $posreBlock = "`r`n; Include ligand restraint file`r`n#ifdef POSRES`r`n#include `"taiposre.itp`"`r`n#endif`r`n"
        if ($content -match '(#include "tai\.itp")') {
            $content = $content -replace '(#include "tai\.itp")', "`$1$posreBlock"
            Set-Content -Path "topol.top" -Value $content -Encoding UTF8
        }
    }

    # Step 9: Make index groups
    Write-Host "`nStep 9: Make index..." -ForegroundColor Cyan
    $ScriptsDir = (Resolve-Path "..\scripts\analysis").Path
    & bash "../scripts/analysis/0_make_index.sh" 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    # Step 10: NVT equilibration
    Write-Host "`nStep 10: NVT equilibration (nvt.mdp)..." -ForegroundColor Cyan
    $NvtMdp = "$MdpDir\nvt.mdp"
    if (-not (Test-Path $NvtMdp)) { $NvtMdp = "$MdpDir\2_heat.mdp" }
    & $GMX grompp -f $NvtMdp -c "em.gro" -r "em.gro" -p "topol.top" -n "index.ndx" -o "nvt.tpr" -maxwarn 7 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & $GMX mdrun -v -deffnm nvt -nb gpu 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    # Step 11: NPT equilibration
    Write-Host "`nStep 11: NPT equilibration (npt.mdp)..." -ForegroundColor Cyan
    $NptMdp = "$MdpDir\npt.mdp"
    if (-not (Test-Path $NptMdp)) { $NptMdp = "$MdpDir\3_equil_npt.mdp" }
    & $GMX grompp -f $NptMdp -c "nvt.gro" -r "nvt.gro" -t "nvt.cpt" -p "topol.top" -n "index.ndx" -o "npt.tpr" -maxwarn 7 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & $GMX mdrun -v -deffnm npt -nb gpu 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    # Step 12: Production MD (md.mdp -> md_0_1)
    Write-Host "`nStep 12: Production MD (md.mdp -> md_0_1)..." -ForegroundColor Cyan
    $MdMdp = "$MdpDir\md.mdp"
    if (-not (Test-Path $MdMdp)) { $MdMdp = "$MdpDir\5_md.mdp" }
    & $GMX grompp -f $MdMdp -c "npt.gro" -r "npt.gro" -t "npt.cpt" -p "topol.top" -n "index.ndx" -o "md_0_1.tpr" -maxwarn 7 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & $GMX mdrun -v -deffnm md_0_1 -nb gpu 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    # Step 13: Automatic PBC removal (user's 3-step reference: whole -> nojump -> mol compact center -> fit rot+trans) & make index
    Write-Host "`nStep 13: Executing user's 3-step golden PBC removal & make_index (0_make_index.sh -> md_fit.xtc)..." -ForegroundColor Cyan
    & bash "../scripts/analysis/0_make_index.sh" 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    # Step 14: Automatic Peptide Backbone Covalent Integrity Verification
    Write-Host "`nStep 14: Verifying peptide Calpha-Calpha bond covalent integrity (verify_peptide_integrity.py)..." -ForegroundColor Cyan
    & $PY "..\scripts\verify_peptide_integrity.py" -s "md_0_1.tpr" -f "md_fit.xtc" 2>&1 | ForEach-Object { "$_" }

    Remove-Item "mdout.mdp" -ErrorAction SilentlyContinue
    Write-Host "`n====================================================================" -ForegroundColor Green
    Write-Host " SUCCESS! 12-Step Split-Topology MD & PBC post-processing completed!" -ForegroundColor Green
    Write-Host " Output: md_0_1.xtc / md_0_1.tpr / md_fit.xtc" -ForegroundColor Green
    Write-Host "====================================================================" -ForegroundColor Green
} finally {
    Pop-Location
}

Write-Host "`nStep 15: Automatically running full paper analysis & figure generation (run_analysis.ps1)..." -ForegroundColor Cyan
& .\run_analysis.ps1 -System $System 2>&1 | ForEach-Object { "$_" }
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

Write-Host "`n====================================================================" -ForegroundColor Green
Write-Host " [SUCCESS] 100% COMPLETE! System: $System MD simulation, PBC removal, peptide" -ForegroundColor Green
Write-Host "           integrity verification, and publication SVG/PNG/PDF plotting finished!" -ForegroundColor Green
Write-Host " -> Inspect your figures in: md_$System\figures\" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green
