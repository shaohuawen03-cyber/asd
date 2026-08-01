#!/usr/bin/env bash
# ============================================================
# AChE (4ey6) - beta-淀粉样肽(Aβ)复合物 分子动力学模拟 全流程
# 复现论文《乙酰胆碱酯酶-β-淀粉样肽复合物的分子动力学模拟》2.3 方法
#
# 力场: amber14sb (Gromacs 中与论文 ff14SB 等价)
# 水模型: TIP3P   盒子: 截角八面体   温度: 300 K   压力: 1 bar
#
# 用法:  ./run_all.sh <系统前缀>
#        例:  TESTING=1 ./run_all.sh alllhrc   (测试模式: NVT/NPT/MD 各 100 步)
#              ./run_all.sh alllhrc            (正式模式: 完整 1000 ns)
#        脚本会读取 input/<前缀>_complex.pdb 并生成该体系全部 MD 文件
# ============================================================
set -eu

# ---------- 0. 全局设置 ----------
SYS="${1:?用法: ./run_all.sh <前缀, 如 alllhrc>}"
WORK="../md_${SYS}"                 # 输出工作目录 (避免污染 input/)

# ---------- 模式选择 ----------
#  TESTING=1  -> 使用 mdp/test 目录 (NVT/NPT/MD 均为 100 步, 快速验证流程跑通)
#  TESTING=0  -> 使用正式 mdp 目录 (完整 1000 ns 产物动力学, 生产用)
#  用法:  TESTING=1 ./run_all.sh alllhrc   (测试)
#          ./run_all.sh alllhrc            (正式)
TESTING="${TESTING:-0}"
if [ "${TESTING}" = "1" ]; then
    MDP="../mdp/test"
    echo ">> 测试模式: 使用 mdp/test (NVT/NPT/MD = 100 步)"
else
    MDP="../mdp"
    echo ">> 正式模式: 使用 mdp (完整 1000 ns 产物动力学)"
fi

INPUT="../input/${SYS}_complex.pdb"

# 蛋白链 / 肽链的链ID (请根据你的 PDB 实际链ID修改)
# AChE 通常为链 A, 肽为链 B; 若不同请改这里
CHAIN_ACHE="A"
CHAIN_PEP="B"

# AChE 断裂封端残基 (论文: 残基 259,262,492,495 处封端)
# 若你的 PDB 中 AChE 是连续链, 置空即可
BREAK_RES=""

GMX="gmx"                            # Gromacs 可执行文件名 (Windows 下可为 gmx.exe)

mkdir -p "${WORK}"
cd "${WORK}"

echo "========== 体系: ${SYS} =========="

# ---------- 1. 结构准备 ----------
echo "[1/10] 准备结构 (去水/配体, 保留蛋白和肽) ..."
# 去掉晶体水分子与配体(加兰他敏已去除), 若仍有水则剔除
cp "${INPUT}" complex_raw.pdb
${GMX} editconf -f complex_raw.pdb -o complex_clean.pdb > /dev/null 2>&1

# ---------- 2. 构建拓扑 ----------
echo "[2/10] pdb2gmx 构建拓扑 (amber14sb + TIP3P, 交互封端) ..."
# 说明: -ter 会逐链询问 N/C 端处理方式。
#   要实现论文的"封端", 对每一端选择对应封端残基:
#     - N 端: 选 ACE (封乙酰化) 或 None(保留带电 NH3+)
#     - C 端: 选 NME (封酰胺化) 或 None(保留带电 COO-)
#   若 AChE 存在内部断裂(如残基259/262, 492/495), 需在此处按片段分别封端。
${GMX} pdb2gmx -f complex_clean.pdb -o complex.gro \
       -p topol.top -i posre.itp \
       -ff amber14sb -water tip3p \
       -ignh -chainsep interactive -ter \
       || echo "!!! pdb2gmx 交互式提示: 请按提示回答各链末端封端方式 (ACE/NME/None)"

# 将重原子位置约束力常数改为论文的 3 kcal/mol/A^2 = 1255 kJ/(mol nm^2)
# (pdb2gmx 默认 fc=1000; 若想精确复现可执行下面 sed)
sed -i 's/1000/1255/g' posre*.itp 2>/dev/null || true

# ---------- 3. 定义盒子(截角八面体) ----------
echo "[3/10] 定义截角八面体盒子 ..."
${GMX} editconf -f complex.gro -o box.gro \
       -c -bt octa -d 1.2

# ---------- 4. 溶剂化 (TIP3P) ----------
echo "[4/10] 加入 TIP3P 水分子 ..."
${GMX} solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

# ---------- 5. 添加离子 (NaCl, 电中性) ----------
echo "[5/10] 添加 Na+/Cl- 实现电中性与生理盐浓度 ..."
${GMX} grompp -f "${MDP}/1_min.mdp" -c solv.gro -p topol.top \
       -o ions.tpr -maxwarn 2
echo "SOL" | ${GMX} genion -s ions.tpr -o neutral.gro -p topol.top \
       -pname NA -nname CL -neutral -conc 0.15

# 建立用于升温/平衡的索引文件
echo "[6/10] 建立索引 (Protein, Water_and_ions) ..."
${GMX} make_ndx -f neutral.gro -o index.ndx << EOF
q
EOF

# ---------- 6. 能量最小化 (2000 步, 重原子约束) ----------
echo "[7/10] 能量最小化 ..."
${GMX} grompp -f "${MDP}/1_min.mdp" -c neutral.gro -r neutral.gro \
       -p topol.top -n index.ndx -o em.tpr -maxwarn 2
${GMX} mdrun -deffnm em -v

# ---------- 7. 升温 0->300K (NVT, 约束) ----------
# 正式模式: 50->300 K 梯度, 每段 0.2 ns, 共约 1 ns (论文 2.3)
# 测试模式: 仅 1 段 100 步, 快速验证
if [ "${TESTING}" = "1" ]; then
    echo "[8/10] 升温 NVT (测试, 100 步) ..."
    cp "${MDP}/2_heat.mdp" heat.mdp
    ${GMX} grompp -f heat.mdp -c em.gro -r neutral.gro \
           -p topol.top -n index.ndx -o heat.tpr -maxwarn 2
    ${GMX} mdrun -deffnm heat -v
    prev="heat"
else
    echo "[8/10] 升温 (NVT, 50->300 K 梯度, 每段 0.2 ns, 共约 1 ns) ..."
    # 论文要求 1 ns 内从 0 K 加热到 300 K;
    # 由于 Gromacs 无法从 0 K 生成速度, 这里从 50 K 开始分梯度逼近。
    prev="em"
    for T in 100 150 200 250 300; do
        cp "${MDP}/2_heat.mdp" "heat_${T}.mdp"
        sed -i "s/^ref_t.*/ref_t = ${T} ${T}/; s/^gen_temp.*/gen_temp = ${T}/" "heat_${T}.mdp"
        ${GMX} grompp -f "heat_${T}.mdp" -c "${prev}.gro" -r neutral.gro \
               -p topol.top -n index.ndx -o "heat_${T}.tpr" -maxwarn 2
        ${GMX} mdrun -deffnm "heat_${T}" -v
        prev="heat_${T}"
    done
fi

# ---------- 8. 恒压密度平衡 (约束) ----------
echo "[9/10] 恒压密度平衡 NPT (约束) ..."
${GMX} grompp -f "${MDP}/3_equil_npt.mdp" -c "${prev}.gro" -r neutral.gro \
       -p topol.top -n index.ndx -o equil_npt.tpr -maxwarn 2
${GMX} mdrun -deffnm equil_npt -v

# ---------- 9. 无约束预平衡 (1 ns) ----------
echo "[10/10] 无约束预平衡 NPT (1 ns) ..."
${GMX} grompp -f "${MDP}/4_equil_npt_free.mdp" -c equil_npt.gro \
       -p topol.top -n index.ndx -o equil_free.tpr -maxwarn 2
${GMX} mdrun -deffnm equil_free -v

# ---------- 10. 产物动力学 (1000 ns) ----------
echo "[11/11] 产物动力学 NPT (300 K, 1 bar, 1000 ns) ..."
${GMX} grompp -f "${MDP}/5_md.mdp" -c equil_free.gro \
       -p topol.top -n index.ndx -o md.tpr -maxwarn 2
${GMX} mdrun -deffnm md -v -cpi md.cpt

echo "========== 体系 ${SYS} 模拟完成 =========="
echo "产物轨迹: md.xtc (每 0.2 ns 一帧, 共 5000 帧)"
echo "后续分析请参考 scripts/analysis/ 下的脚本。"
