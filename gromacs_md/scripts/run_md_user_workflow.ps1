# ============================================================
# GROMACS 100 ns MD Workflow following User's Standard Protocol & mdp_templates
#
# Usage:  .\run_md_user_workflow.ps1 -System alllhrc
#         .\run_md_user_workflow.ps1 -System fllhttr
#         .\run_md_user_workflow.ps1 -System ylsllqr
#         .\run_md_user_workflow.ps1 -System ache
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

# Locate mdp_templates directory (checking standard user path F:\0wsh\asd\mdp_templates or relative paths)
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

# Resolve PDB input file
$InputPdb = ""
foreach ($pdbCand in @("..\input\${System}_complex.pdb", "..\input\${System}.pdb", "..\input\${System}_alone.pdb")) {
    if (Test-Path $pdbCand) {
        $InputPdb = (Resolve-Path $pdbCand).Path
        break
    }
}

# If standalone ache monomer requested and not present, automatically extract from alllhrc_complex.pdb
if (-not $InputPdb -and ($System -eq "ache" -or $System -eq "ache_alone" -or $System -eq "pro")) {
    Write-Host ">> Automatically extracting standalone AChE monomer from alllhrc_complex.pdb..." -ForegroundColor Yellow
    & $PY ".\extract_ache_monomer.py" "..\input\alllhrc_complex.pdb" "..\input\ache.pdb" 2>&1 | ForEach-Object { "$_" }
    if (Test-Path "..\input\ache.pdb") {
        $InputPdb = (Resolve-Path "..\input\ache.pdb").Path
    }
}
if (-not $InputPdb) {
    Write-Host "ERROR: Could not find input PDB file for System: $System!" -ForegroundColor Red
    exit 1
}

Write-Host "====================================================================" -ForegroundColor Green
Write-Host " Starting 100 ns MD Simulation for System: $System (User Workflow)" -ForegroundColor Green
Write-Host " Input PDB: $InputPdb" -ForegroundColor Green
Write-Host "====================================================================" -ForegroundColor Green

Push-Location $WorkDir
try {
    # Step 1: Generate topology (pdb2gmx)
    Write-Host "`nStep 1: Running gmx pdb2gmx..." -ForegroundColor Cyan
    Copy-Item $InputPdb "pro.pdb" -Force
    @("5", "5") | & $GMX pdb2gmx -f "pro.pdb" -o "pro_processed.gro" -p "topol.top" -ignh -ff "amber99sb-ildn" -water "tip3p" 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    # Step 4: editconf - define triclinic box (-d 1.0)
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

    # Step 12: Production MD (100 ns -> md_0_1)
    Write-Host "`nStep 12: Production MD (md.mdp -> md_0_1)..." -ForegroundColor Cyan
    $MdMdp = "$MdpDir\md.mdp"
    if (-not (Test-Path $MdMdp)) { $MdMdp = "$MdpDir\5_md.mdp" }
    & $GMX grompp -f $MdMdp -c "npt.gro" -r "npt.gro" -t "npt.cpt" -p "topol.top" -n "index.ndx" -o "md_0_1.tpr" -maxwarn 7 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & $GMX mdrun -v -deffnm md_0_1 -nb gpu 2>&1 | ForEach-Object { "$_" }
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    # Clean intermediate mdout.mdp
    Remove-Item "mdout.mdp" -ErrorAction SilentlyContinue

    Write-Host ""
    Write-Host "====================================================================" -ForegroundColor Green
    Write-Host " SUCCESS! System $System 100 ns MD completed! Output: md_0_1.xtc / md_0_1.tpr" -ForegroundColor Green
    Write-Host "====================================================================" -ForegroundColor Green
} finally {
    Pop-Location
}
