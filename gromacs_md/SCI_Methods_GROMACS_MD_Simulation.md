# SCI 论文“材料与方法”完整规范描述：乙酰胆碱酯酶–Aβ 肽对接复合物 100 ns 分子动力学模拟与多层次规范分析

**Materials and Methods: Molecular Dynamics Simulations and Multi-Level Trajectory Analyses of AChE-Aβ Complexes**

---

## 1. 起始体系准备与分子对接结构处理 (System Preparation and Docking Conformations)
为了探究不同肽段（如 `alllhrc`、`fllhttr`、`ylsllqr` 7肽构象及参考 `Aβ(1-42)` 肽）与乙酰胆碱酯酶（Acetylcholinesterase, AChE）外周阴离子位点（Peripheral Anionic Site, PAS）的相互作用机制，所有起始复合物结构基于分子对接的最优构象构建。在进行分子动力学（MD）模拟前，对于目标 AChE（参考晶体结构 PDB ID: **4ey6**，分辨率 2.40 Å），依据论文方法要求，对蛋白内部断裂部位（残基 259/262、492/495）以及链 N端 / C端执行末端处理（如 `-ter` 标志交互分配乙酰化 `ACE` 封端或 N-甲基酰胺化 `NME` 封端，或按照默认生成标准带电末端 `NH3+` / `COO-`）。所有结晶水分子及杂质分子被事先剥离以重建均一的显式溶剂环境。

## 2. 力场选择、三斜盒子构建与生理盐溶剂化 (Force Field, Triclinic Box & Solvation)
所有分子动力学模拟基于 **GROMACS 2025** 软件套件在等压等温（NPT）及等容等温（NVT）系综下完成。力场采用与 **AMBER ff14SB** 完全等价的标准推荐力场 **`amber99sb-ildn`**（或者可选手动导入的 `amber14sb.ff` 社区移植版）。

> **★ 三斜盒子优化说明**：为了最大限度地节约溶剂化体积与非键相互作用运算耗时，同时严格保证周期性边界条件（PBC）下最小映射镜像距离大于 2.0 nm，采用三斜周期性盒子（Triclinic Periodic Box, `gmx editconf -c -bt triclinic -d 1.0`），设定溶质任一表面原子距离三斜盒边界的几何垂直间距为 1.0 nm (`-d 1.0`)。与传统截角八面体（Octahedral, `d=1.2 nm`）相比，三斜盒子显著缩减了 **~35%** 的冗余溶剂化水分子体积。

系统随后在三斜盒子中加入标准 **TIP3P** 显式水模型 (`spc216.gro`)，并通过 `gmx genion` 以 Joung-Cheatham 离子参数体系添加足量的对等离子与适量的 `Na+` 和 `Cl-`，中和系统总电荷并确保体系中生理盐浓度达到 **0.15 mol/L (0.15 M NaCl)**。

## 3. 能量最小化与多阶段约束平衡模拟 (EM and Multi-Stage Equilibration Protocol)
在正式产物模拟开始前，对溶剂化中和体系执行严格的阶梯式弛豫与平衡步骤：
1. **能量最小化（EM）**：使用最速下降法 (`steep`) 在对所有蛋白质和肽的重原子施加位置约束（`define = -DPOSRES`，力常数 $k = 3 \text{ kcal} \cdot \text{mol}^{-1} \cdot \text{\AA}^{-2} = 1255 \text{ kJ} \cdot \text{mol}^{-1} \cdot \text{nm}^{-2}$）的前提下进行 2,000 步能量优化，消除位阻冲突。
2. **NVT 升温模拟（0 K 至 300 K，1.0 ns）**：将 1.0 ns 的 NVT 加热过程划分成 5 个连续的等时分段（100 K -> 150 K -> 200 K -> 250 K -> 300 K，每个温段 0.2 ns）。升温阶段对所有重原子保留 1255 $\text{kJ}/(\text{mol} \cdot \text{nm}^2)$ 位置约束，控温采用 Velocity-rescale (`v-rescale`)，时间常数 $\tau_T = 1.0 \text{ ps}$。
3. **NPT 恒压密度平衡（1.0 ns，重原子约束）**：继续保留重原子约束，在 300 K 和 1.0 bar 压力下运行 1.0 ns (`500,000` 步)。控压采用 Berendsen 压力耦合器（`isotropic`，$\tau_P = 2.0 \text{ ps}$，压缩率 $4.5 \times 10^{-5} \text{ bar}^{-1}$）。
4. **NPT 无约束预平衡（1.0 ns）**：释放全体系所有位置约束（`define = -DFLEXIBLE`），在 300 K、1.0 bar 条件下继续运行 1.0 ns 无约束预平衡。

## 4. 100 ns 产物分子动力学模拟 (100 ns Production MD Simulation)
各复合物体系在 NPT 系综（300 K, 1.0 bar）下完成 **100 ns** 的正式产物分子动力学模拟（Production MD）。步长 $\text{dt} = 2.0 \text{ fs}$，100 ns 共执行 **`50,000,000` 步**。模拟积分采用 Leap-frog 积分算法。对涉及氢原子的共价键采用 **LINCS** 算法约束。非键相互作用采用 **Verlet** 截断机制，截断半径 $r = 1.2 \text{ nm}$（范德华力平滑转换从 1.0 nm 开始）。长程静电相互作用通过 **PME (Particle Mesh Ewald)** 方法计算。系统每隔 **20 ps (`0.02 ns`, `nstxout-compressed = 10000`)** 记录一个坐标帧， entire 100 ns 轨迹共收集 **5,000 个分析构象帧**。

## 5. 模拟轨迹对标论文多层次分析方法 (Trajectory Analysis Methodology)
对 100 ns 模拟轨迹执行的完整分析工作流对应于论文的 **图 1–图 6 与 表 1–表 2**：
- **骨架 RMSD 与 RMSF（图 1，论文 3.1 节）**：使用 `gmx rms` 和 `gmx rmsf` 计算复合物整体、AChE 与肽链骨架 Cα 的 RMSD 与 RMSF。
- **质心径向分布函数 RDF（图 2，论文 3.1 节）**：使用 `gmx rdf -selrpos mol_com -seltype mol_com` 计算肽分子围绕 AChE 分子质心的分布，分 4 等份验证收敛性。
- **溶剂可及表面积 SASA（图 3，论文 3.2 节）**：基于类 Shrake-Rupley 算法（LCPO 等效算法，`gmx sasa`）计算 SASA 均值。
- **二级结构动态分布（图 4，论文 3.2 节）**：通过 GROMACS `dssp` 模块逐帧识别并计算每 50 ns 内 α-螺旋、转角、弯曲的分布。
- **分子间与内部氢键统计（论文 3.3 节）**：通过 `gmx hbond`（距离 < 3.0 Å，角度 30° / 135°）记录氢键时间序列。
- **天然与非天然残基对接触统计（图 5 / 表 1，论文 3.3 节）**：基于 Python MDAnalysis 计算残基间距矩阵（< 7.0 Å 判据），识别频次超过 10 次的非天然残基对。
- **水介导桥连相互作用统计（图 6 / 表 2，论文 3.4 节）**：统计同时处于肽和 AChE 势能范围内的显式桥连水分子 (`SOL`)，记录逐残基桥连水数与作用数。

---

## 6. 方法学配置对照表 (Reference Tables)

### 表 1. 100 ns 动力学模拟全阶段 MDP 参数对照表
| 模拟阶段 | 时长/步数 | 积分步长 (dt) | 约束条件 | 控温 / 控压方式 |
|---|---|---|---|---|
| **1. 能量最小化 (EM)** | 2,000 步 (Fmax < 1000) | 0.01 (emstep) | 重原子 1255 kJ/mol/nm² | 最速下降法 (Steepest Descents) |
| **2. NVT 连续退火升温** | 1.0 ns (500,000步) | 2.0 fs (0.002) | 重原子 1255 kJ/mol/nm² | 模拟退火 (0->300 K 连续线性升温) |
| **3. NPT 恒压平衡** | 1.0 ns (500,000步) | 2.0 fs (0.002) | 重原子 1255 kJ/mol/nm² | v-rescale (300K) / Berendsen (tau_P=5.0 ps) |
| **4. NPT 无约束预平衡** | 1.0 ns (500,000步) | 2.0 fs (0.002) | 无 (-DFLEXIBLE) | v-rescale (300K) / Berendsen (tau_P=5.0 ps) |
| **5. 正式产物模拟 (MD)** | 100 ns (50,000,000步) | 2.0 fs (0.002) | 无约束 / LINCS 氢键 | v-rescale (300K) / Berendsen/PR (tau_P=5.0 ps) |

### 表 2. 模拟产物轨迹分析命令与论文图表严格映射对照表
| 分析项目 / 指标 | 论文章节与对应图表 | 核心执行指令或算法脚本 | 输出产物文件 |
|---|---|---|---|
| **三斜盒子定义** | 论文 2.2 / 节约计算体积 | `gmx editconf -c -bt triclinic -d 1.0` | `box.gro` (体积减少 ~35%) |
| **骨架 RMSD / RMSF** | 论文 3.1 节 / 图 1A–1F | `gmx rms` & `gmx rmsf (-fit rot+trans)` | `rmsd_*_bb.xvg` / `rmsf_*_bb.xvg` |
| **径向分布函数 (RDF)** | 论文 3.1 节 / 图 2A–2B | `gmx rdf -selrpos mol_com -seltype mol_com` | `rdf_pep_ache.xvg` / `Q1-Q4.xvg` |
| **溶剂可及表面积 SASA** | 论文 3.2 节 / 图 3A–3B | `gmx sasa -surface Protein -output Protein` | `sasa_complex.xvg` / `sasa_pep.xvg` |
| **二级结构演变 (DSSP)** | 论文 3.2 节 / 图 4 | `gmx dssp -sel Peptide` & `dssp_bins.py` | `ss_pep.dat` / `ss_pep_bins.dat` |
| **氢键数量动态曲线** | 论文 3.3 节 | `gmx hbond` (几何条件 < 3.0 Å, 角 < 30°) | `hbond_ache_pep.xvg` |
| **非天然残基对接触** | 论文 3.3 节 / 图 5 & 表 1 | `contacts.py (-cut 7.0 Å, MDAnalysis)` | `inter_contacts.csv` / `frequent_contacts.tsv` |
| **水介导桥连相互作用** | 论文 3.4 节 / 图 6 & 表 2 | `bridging_waters.py` (近邻水过滤 100x 加速) | `bridging_per_residue.csv` |
| **全套一键批量图表生成** | 全面对应所有发表级图表 | `plot_all.py` (自动添加 A, B, C... 子图编号) | `fig0_summary_all.svg/.png/.pdf` |
