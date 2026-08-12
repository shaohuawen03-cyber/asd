#!/usr/bin/env bash
# ============================================================
# GROMACS 100 ns MD Workflow following User's Standard Protocol & mdp_templates (Bash / WSL)
#
# Usage:  ./run_md_user_workflow.sh alllhrc
#         ./run_md_user_workflow.sh fllhttr
#         ./run_md_user_workflow.sh ylsllqr
#         ./run_md_user_workflow.sh ache
# ============================================================
set -eu

SYS="${1:?用法: ./run_md_user_workflow.sh <系统前缀, 如 alllhrc>}"

GMX="gmx"
if command -v gmx.exe >/dev/null 2>&1; then
    GMX="gmx.exe"
fi

PY="python"
if command -v python3 >/dev/null 2>&1; then
    PY="python3"
fi

# Locate mdp_templates directory
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
echo " Starting 100 ns MD Simulation for System: ${SYS} (User Workflow)"
echo " Input PDB: ${INPUT_PDB}"
echo "===================================================================="

cd "${WORK}"

echo ""
echo "Step 1: Running gmx pdb2gmx..."
cp "${INPUT_PDB}" "pro.pdb"
echo -e "5\n5" | ${GMX} pdb2gmx -f "pro.pdb" -o "pro_processed.gro" -p "topol.top" -ignh -ff "amber99sb-ildn" -water "tip3p"

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

rm -f "mdout.mdp" 2>/dev/null || true

echo ""
echo "===================================================================="
echo " SUCCESS! System ${SYS} 100 ns MD completed! Output: md_0_1.xtc / md_0_1.tpr"
echo "===================================================================="
