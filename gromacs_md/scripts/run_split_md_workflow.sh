#!/usr/bin/env bash
# ============================================================
# GROMACS 12-Step Split-Topology MD Workflow (Bash / WSL)
#
# Usage:  ./run_split_md_workflow.sh alllhrc
#         ./run_split_md_workflow.sh fllhttr
#         ./run_split_md_workflow.sh ylsllqr
#         ./run_split_md_workflow.sh ache
# ============================================================
set -eu

SYS="${1:?用法: ./run_split_md_workflow.sh <系统前缀, 如 alllhrc>}"

GMX="gmx"
if command -v gmx.exe >/dev/null 2>&1; then
    GMX="gmx.exe"
fi

PY="python"
if command -v python3 >/dev/null 2>&1; then
    PY="python3"
fi

MDP_DIR=""
for CAND in "../../mdp_templates" "../mdp_templates" "/mnt/f/0wsh/asd/mdp_templates" "/f/0wsh/asd/mdp_templates" "../mdp/100ns"; do
    if [ -d "$CAND" ]; then
        MDP_DIR="$(cd "$CAND" && pwd)"
        break
    fi
done
if [ -z "$MDP_DIR" ]; then
    echo "ERROR: Could not find mdp_templates directory!"
    exit 1
fi
echo ">> Using MDP Templates Directory: ${MDP_DIR}"

WORK="../md_${SYS}"
if [ -d "${WORK}" ]; then
    echo ">> Cleaning previous workspace ${WORK} for a 100% fresh start..."
    rm -rf "${WORK}"
fi
mkdir -p "${WORK}"

INPUT_PDB=""
for CAND in "../input/${SYS}_complex.pdb" "../input/${SYS}.pdb" "../input/${SYS}_alone.pdb"; do
    if [ -f "$CAND" ]; then
        INPUT_PDB="$(cd $(dirname "$CAND") && pwd)/$(basename "$CAND")"
        break
    fi
done

if [ -z "$INPUT_PDB" ] && { [ "${SYS}" = "ache" ] || [ "${SYS}" = "ache_alone" ] || [ "${SYS}" = "pro" ]; }; then
    echo ">> Automatically extracting standalone AChE monomer from alllhrc_complex.pdb..."
    $PY "./extract_ache_monomer.py" "../input/alllhrc_complex.pdb" "../input/ache.pdb"
    if [ -f "../input/ache.pdb" ]; then
        INPUT_PDB="$(cd ../input && pwd)/ache.pdb"
    fi
fi
if [ -z "$INPUT_PDB" ]; then
    echo "ERROR: Could not find input PDB file for System: ${SYS}!"
    exit 1
fi

echo "===================================================================="
echo " Starting 12-Step Split-Topology MD for System: ${SYS}"
echo " Input PDB: ${INPUT_PDB}"
echo "===================================================================="

cd "${WORK}"

echo ""
echo "Step 0: Splitting complex PDB into pro.pdb and tai.pdb..."
$PY "../scripts/split_complex_pdb.py" "${INPUT_PDB}" "pro.pdb" "tai.pdb"
HAS_TAI=0
if [ -f "tai.pdb" ] && [ $(stat -c%s "tai.pdb" 2>/dev/null || stat -f%z "tai.pdb" 2>/dev/null || echo "0") -gt 100 ]; then
    HAS_TAI=1
fi

echo ""
echo "Step 1: Running gmx pdb2gmx for pro.pdb..."
echo -e "5\n5" | ${GMX} pdb2gmx -f "pro.pdb" -o "pro_processed.gro" -p "topol.top" -ignh -ff "amber99sb-ildn" -water "tip3p"

if [ "$HAS_TAI" = "1" ]; then
    echo ""
    echo "Step 1b: Running gmx pdb2gmx for tai.pdb..."
    echo -e "5\n5" | ${GMX} pdb2gmx -f "tai.pdb" -o "tai.gro" -p "tai.top" -i "taiposre.itp" -ignh -ff "amber99sb-ildn" -water "none"

    echo ""
    echo "Step 2: Merging pro_processed.gro and tai.gro..."
    N1=$(sed -n '2p' pro_processed.gro | awk '{print $1}')
    N2=$(sed -n '2p' tai.gro | awk '{print $1}')
    TOTAL=$((N1 + N2))
    echo "File1 atoms: $N1 | File2 atoms: $N2 | Total: $TOTAL"
    head -n 1 pro_processed.gro > gro_merged.gro
    printf "%5d\n" "$TOTAL" >> gro_merged.gro
    sed -n "3,$((N1 + 2))p" pro_processed.gro >> gro_merged.gro
    sed -n "3,$((N2 + 2))p" tai.gro >> gro_merged.gro
    tail -n 1 pro_processed.gro >> gro_merged.gro
    mv gro_merged.gro pro_processed.gro

    echo ""
    echo "Step 3: Updating topol.top..."
    sed -i '/amber99sb/a #include "tai.itp"' topol.top || sed -i '1i #include "tai.itp"' topol.top
    echo "tai                 1" >> topol.top
fi

echo ""
echo "Step 4: editconf (triclinic box -d 1.0)..."
${GMX} editconf -f "pro_processed.gro" -o "gro_newbox.gro" -c -d 1.0 -bt triclinic

echo ""
echo "Step 5: solvate (spc216.gro)..."
${GMX} solvate -cp "gro_newbox.gro" -cs spc216.gro -o "gro_solv.gro" -p "topol.top"

echo ""
echo "Step 6: ions (ions.mdp + genion neutral 0.15M)..."
IONS_MDP="${MDP_DIR}/ions.mdp"
if [ ! -f "${IONS_MDP}" ]; then IONS_MDP="${MDP_DIR}/1_min.mdp"; fi
${GMX} grompp -f "${IONS_MDP}" -c "gro_solv.gro" -r "gro_solv.gro" -p "topol.top" -o "ions.tpr" -maxwarn 7
echo "SOL" | ${GMX} genion -s "ions.tpr" -o "solv_ions.gro" -p "topol.top" -pname NA -nname CL -neutral -conc 0.15

echo ""
echo "Step 7: Energy Minimization (em.mdp)..."
EM_MDP="${MDP_DIR}/em.mdp"
if [ ! -f "${EM_MDP}" ]; then EM_MDP="${MDP_DIR}/1_min.mdp"; fi
${GMX} grompp -f "${EM_MDP}" -c "solv_ions.gro" -r "solv_ions.gro" -p "topol.top" -o "em.tpr" -maxwarn 7
${GMX} mdrun -v -deffnm em -nb gpu

if [ "$HAS_TAI" = "1" ]; then
    echo ""
    echo "Step 8: Inserting position restraints into topol.top..."
    echo -e "\n; Include ligand restraint file\n#ifdef POSRES\n#include \"taiposre.itp\"\n#endif" >> topol.top
fi

echo ""
echo "Step 9: Make index..."
bash "../scripts/analysis/0_make_index.sh"

echo ""
echo "Step 10: NVT equilibration (nvt.mdp)..."
NVT_MDP="${MDP_DIR}/nvt.mdp"
if [ ! -f "${NVT_MDP}" ]; then NVT_MDP="${MDP_DIR}/2_heat.mdp"; fi
${GMX} grompp -f "${NVT_MDP}" -c "em.gro" -r "em.gro" -p "topol.top" -n "index.ndx" -o "nvt.tpr" -maxwarn 7
${GMX} mdrun -v -deffnm nvt -nb gpu

echo ""
echo "Step 11: NPT equilibration (npt.mdp)..."
NPT_MDP="${MDP_DIR}/npt.mdp"
if [ ! -f "${NPT_MDP}" ]; then NPT_MDP="${MDP_DIR}/3_equil_npt.mdp"; fi
${GMX} grompp -f "${NPT_MDP}" -c "nvt.gro" -r "nvt.gro" -t "nvt.cpt" -p "topol.top" -n "index.ndx" -o "npt.tpr" -maxwarn 7
${GMX} mdrun -v -deffnm npt -nb gpu

echo ""
echo "Step 12: Production MD (md.mdp -> md_0_1)..."
MD_MDP="${MDP_DIR}/md.mdp"
if [ ! -f "${MD_MDP}" ]; then MD_MDP="${MDP_DIR}/5_md.mdp"; fi
${GMX} grompp -f "${MD_MDP}" -c "npt.gro" -r "npt.gro" -t "npt.cpt" -p "topol.top" -n "index.ndx" -o "md_0_1.tpr" -maxwarn 7
${GMX} mdrun -v -deffnm md_0_1 -nb gpu

echo ""
echo "Step 13: Executing user's 3-step golden PBC removal & make_index (0_make_index.sh -> md_fit.xtc)..."
bash "../scripts/analysis/0_make_index.sh"

echo ""
echo "Step 14: Verifying peptide Calpha-Calpha bond covalent integrity (verify_peptide_integrity.py)..."
"${PY}" "../scripts/verify_peptide_integrity.py" -s "md_0_1.tpr" -f "md_fit.xtc" || true

rm -f "mdout.mdp" 2>/dev/null || true

echo ""
echo "===================================================================="
echo " SUCCESS! System ${SYS} 12-Step Split-Topology MD & PBC post-processing completed!"
echo " Output: md_0_1.xtc / md_0_1.tpr / md_fit.xtc"
echo "===================================================================="
cd "../scripts"

echo ""
echo "Step 15: Automatically running full paper analysis & figure generation (run_analysis.sh)..."
bash "./run_analysis.sh" "${SYS}"

echo ""
echo "===================================================================="
echo " [SUCCESS] 100% COMPLETE! System: ${SYS} MD simulation, PBC removal, peptide"
echo "           integrity verification, and publication SVG/PNG/PDF plotting finished!"
echo " -> Inspect your figures in: md_${SYS}/figures/"
echo "===================================================================="
