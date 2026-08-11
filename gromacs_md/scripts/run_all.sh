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

# ---------- 模拟参数配置 ----------
#   SIM_MODE="100ns"   -> 100 ns 正式产物动力学模拟 (默认配置: mdp/100ns)
#   SIM_MODE="1000ns"  -> 1000 ns 正式产物动力学模拟 (mdp)
SIM_MODE="${SIM_MODE:-100ns}"
if [ "${SIM_MODE}" = "100ns" ] || [ "${TESTING}" = "1" ]; then
    MDP="../mdp/100ns"
    echo ">> [100 ns 正式生产模拟] 使用 mdp/100ns (NVT梯度 1.0 ns -> NPT 2.0 ns -> MD 100 ns -> 5,000 帧)"
else
    MDP="../mdp"
    echo ">> [1000 ns 正式生产模拟] 使用 mdp (完整 1000 ns 产物动力学)"
fi

INPUT="../input/${SYS}_complex.pdb"

# 蛋白链 / 肽链的链ID (请根据你的 PDB 实际链ID修改)
# AChE 通常为链 A, 肽为链 B; 若不同请改这里
CHAIN_ACHE="A"
CHAIN_PEP="B"

# AChE 断裂封端残基 (论文: 残基 259,262,492,495 处封端)
# 若你的 PDB 中 AChE 是连续链, 置空即可
BREAK_RES=""

# ---------- 开放配置参数 ----------
# 1. 力场选择 (标准 GROMACS 默认自带 amber99sb-ildn;
#   若你在本地安装了第三方 amber14sb.ff 力场包, 可在命令中加入 FORCE_FIELD=amber14sb)
FORCE_FIELD="${FORCE_FIELD:-amber99sb-ildn}"

# 2. 链末端封端交互设置 (-ter)
#   INTERACTIVE_TER=0 -> 全自动无交互，自动使用标准的带电 N/C 末端 (默认, 适合自动化跑通)
#   INTERACTIVE_TER=1 -> 开启 -ter 交互模式，在命令行询问 N/C 端封端方式 (ACE/NME/None)
INTERACTIVE_TER="${INTERACTIVE_TER:-0}"
TER_FLAG=""
if [ "${INTERACTIVE_TER}" = "1" ]; then
    TER_FLAG="-ter"
fi

# 3. GPU 加速参数 (默认开启 -nb gpu; 若需自定义可传入 GPU_FLAGS="-nb gpu -pme gpu")
GPU_FLAGS="${GPU_FLAGS:--nb gpu}"

# ---------- 自动检测/校验 GROMACS 命令 ----------
if command -v gmx.exe >/dev/null 2>&1; then
    GMX="gmx.exe"
elif command -v gmx >/dev/null 2>&1; then
    GMX="gmx"
elif [ -n "${GROMACS_CMD:-}" ]; then
    GMX="${GROMACS_CMD}"
else
    echo "!!! 错误: 在当前的 Bash 环境的 PATH 环境变量中找不到 'gmx' 或 'gmx.exe' 命令！"
    echo "!!! 常见原因:"
    echo "!!!   1. GROMACS 装在 Windows 目录下（如 C:\\Program Files\\Gromacs\\bin）但没有加入到 PATH 中。"
    echo "!!!   2. 你使用的是 Conda 环境中的 GROMACS，但通过 PowerShell 执行 bash -c 时没有带入 PATH。"
    echo "!!! 解决方案:"
    echo "!!!   - 将 gmx.exe 所在目录放入环境变量 PATH 中；"
    echo "!!!   - 或直接指定命令完整路径: GROMACS_CMD='/c/Program Files/Gromacs/bin/gmx.exe' TESTING=1 ./run_all.sh alllhrc"
    exit 1
fi
echo ">> 使用 GROMACS 执行命令: $(command -v "${GMX}" 2>/dev/null || echo "${GMX}")"

mkdir -p "${WORK}"
cd "${WORK}"

echo "========== 体系: ${SYS} =========="

# ----- 支持单独从第11步产物MD启动 (节省时间, 无需重跑前期 10 步平衡) -----
if [ "${ONLY_MD:-0}" = "1" ] || [ "${2:-}" = "--only-md" ] || [ "${2:-}" = "-OnlyMD" ] || [ "${2:-}" = "onlymd" ]; then
    echo ">> [仅启动产物模拟模式] 跳过前 10 步准备与平衡，直接启动第 [11/11] 步产物 MD..."
    if [ ! -f "equil_free.gro" ]; then
        echo "!!! 错误: 未在工作目录找到 equil_free.gro！请先运行完整的预平衡。"
        exit 1
    fi
    # 删掉原子数不匹配的旧残留 md.cpt，保证能从新盒子纯净启动
    rm -f md.cpt md.part*.cpt 2>/dev/null || true
    ${GMX} grompp -f "${MDP}/5_md.mdp" -c equil_free.gro -r equil_free.gro \
           -p topol.top -n index.ndx -o md.tpr -maxwarn 2
    ${GMX} mdrun -deffnm md -v ${GPU_FLAGS}
    echo "========== 体系 ${SYS} 正式产物模拟完成 =========="
    exit 0
fi

# ---------- 1. 结构准备 ----------
echo "[1/10] 准备结构 (检查输入文件) ..."
if [ ! -f "${INPUT}" ]; then
    echo "!!! 错误: 找不到输入 PDB 文件 '${INPUT}'"
    echo "!!! 请确认你已将 ${SYS}_complex.pdb 放在了 gromacs_md/input/ 目录下"
    exit 1
fi
cp "${INPUT}" complex_clean.pdb

# ---------- 2. 构建拓扑 ----------
echo "[2/10] pdb2gmx 构建拓扑 (力场: ${FORCE_FIELD}, 水模型: TIP3P) ..."
# 说明: -ter 会逐链询问 N/C 端处理方式。
#   要实现论文的"封端", 对每一端选择对应封端残基:
#     - N 端: 选 ACE (封乙酰化) 或 None(保留带电 NH3+)
#     - C 端: 选 NME (封酰胺化) 或 None(保留带电 COO-)
#   若 AChE 存在内部断裂(如残基259/262, 492/495), 需在此处按片段分别封端。
if [ "${INTERACTIVE_TER}" = "1" ]; then
    echo ">> 提示: 已开启交互式封端(-ter)，请在屏幕提示时输入数字选择每一条链的 N/C 端封端方式(如 ACE/NME/None)。"
else
    echo ">> 提示: 自动末端处理模式(未开启-ter)，将默认采用带电 N/C 端。如需交互选择可传入 INTERACTIVE_TER=1。"
fi
${GMX} pdb2gmx -f complex_clean.pdb -o complex.gro \
       -p topol.top -i posre.itp \
       -ff "${FORCE_FIELD}" -water tip3p \
       -ignh ${TER_FLAG}

# 将重原子位置约束力常数改为论文的 3 kcal/mol/A^2 = 1255 kJ/(mol nm^2)
# (pdb2gmx 默认 fc=1000; 若想精确复现可执行下面 sed)
sed -i 's/1000/1255/g' posre*.itp 2>/dev/null || true

# ---------- 3. 定义三斜盒子 (triclinic, d = 1.0 nm) ----------
echo "[3/10] 定义三斜盒子 (triclinic, -d 1.0 nm 节省溶剂化体积与计算耗时) ..."
${GMX} editconf -f complex.gro -o box.gro \
       -c -bt triclinic -d 1.0

# ---------- 4. 溶剂化 (TIP3P) ----------
echo "[4/10] 加入 TIP3P 水分子 ..."
${GMX} solvate -cp box.gro -cs spc216.gro -o solv.gro -p topol.top

# ---------- 5. 添加离子 (NaCl, 电中性) ----------
echo "[5/10] 添加 Na+/Cl- 实现电中性与生理盐浓度 ..."
${GMX} grompp -f "${MDP}/1_min.mdp" -c solv.gro -r solv.gro -p topol.top \
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
${GMX} mdrun -deffnm em -v ${GPU_FLAGS}

# ---------- 7. 升温 0->300K (NVT, 约束, 连续线性模拟退火共 1.0 ns) ----------
echo "[8/10] 升温 (NVT, 0->300 K 连续模拟退火升温 1.0 ns, 严格复现论文) ..."
${GMX} grompp -f "${MDP}/2_heat.mdp" -c em.gro -r neutral.gro \
       -p topol.top -n index.ndx -o heat.tpr -maxwarn 2
${GMX} mdrun -deffnm heat -v ${GPU_FLAGS}
prev="heat"

# ---------- 8. 恒压密度平衡 (约束) ----------
echo "[9/10] 恒压密度平衡 NPT (约束) ..."
${GMX} grompp -f "${MDP}/3_equil_npt.mdp" -c "${prev}.gro" -r neutral.gro \
       -p topol.top -n index.ndx -o equil_npt.tpr -maxwarn 2
${GMX} mdrun -deffnm equil_npt -v ${GPU_FLAGS}

# ---------- 9. 无约束预平衡 (1 ns) ----------
echo "[10/10] 无约束预平衡 NPT (1 ns) ..."
${GMX} grompp -f "${MDP}/4_equil_npt_free.mdp" -c equil_npt.gro -r equil_npt.gro \
       -p topol.top -n index.ndx -o equil_free.tpr -maxwarn 2
${GMX} mdrun -deffnm equil_free -v ${GPU_FLAGS}

# ---------- 10. 产物动力学 (100 ns / 1000 ns) ----------
echo "[11/11] 产物动力学 NPT (300 K, 1 bar, 正式产物动力学) ..."
${GMX} grompp -f "${MDP}/5_md.mdp" -c equil_free.gro -r equil_free.gro \
       -p topol.top -n index.ndx -o md.tpr -maxwarn 2
if [ -f "md.cpt" ]; then
    echo ">> 尝试从 existing md.cpt 续跑产物动力学..."
    ${GMX} mdrun -deffnm md -v -cpi md.cpt ${GPU_FLAGS} || {
        echo ">> [提示] 检测到原有 md.cpt 与当前新体系原子数/三斜盒子不匹配，自动清理旧 cpt 纯净启动 MD..."
        rm -f md.cpt md.part*.cpt
        ${GMX} mdrun -deffnm md -v ${GPU_FLAGS}
    }
else
    ${GMX} mdrun -deffnm md -v ${GPU_FLAGS}
fi

echo "========== 体系 ${SYS} 模拟完成 =========="
echo "产物轨迹: md.xtc (每 0.2 ns 一帧, 共 5000 帧)"
echo "后续分析请参考 scripts/analysis/ 下的脚本。"
