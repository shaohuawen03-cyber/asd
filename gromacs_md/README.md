# AChE(4ey6)–β-淀粉样肽 复合物 分子动力学模拟流程

复现论文《乙酰胆碱酯酶-β-淀粉样肽复合物的分子动力学模拟》的完整 **Gromacs** 工作流。
针对你的三个对接复合物（`alllhrc`、`fllhttr`、`ylsllqr`，均为不同肽与 **4ey6** 在 **PAS 位点** 的最优对接构象）生成测试动力学代码。

## 0. 依赖与准备

- **Gromacs** ≥ 2018（默认使用官方自带的 `amber99sb-ildn` 力场与 `spc216.gro`；如需第三方 `amber14sb.ff` 请传入 `FORCE_FIELD=amber14sb`）
- 分析脚本需要 **Python 3** + **MDAnalysis** + **numpy**：`pip install numpy MDAnalysis`
- 将你的三个 PDB 放入 `input/`，命名方式：
  - `input/alllhrc_complex.pdb`
  - `input/fllhttr_complex.pdb`
  - `input/ylsllqr_complex.pdb`

> Windows 提示：Gromacs 在 Windows 下建议用 WSL/原生 Linux 运行；本套脚本为 bash 脚本，
> 在 WSL 中执行 `./scripts/run_all.sh <前缀>` 即可。`GMX` 变量在 `run_all.sh` 顶部，Windows 下可改为 `gmx.exe`。

## 1. 目录结构

```
gromacs_md/
├── input/                       # 放置三个 *_complex.pdb
├── mdp/                         # ★ 正式(1000 ns 生产)参数
│   ├── 100ns/                   # ★ 100 ns 正式模拟参数配置 (采用 0->300K 连续线性退火)
│   ├── test/                    # ★ 测试参数 (已更新同步为 100 ns 验证套件)
├── scripts/
│   ├── run_100ns_formal.ps1/.sh # ★ 100 ns 正式产物模拟一键全自动套件 (自动清除测试残留 -> 100ns MD -> 绘图)
│   ├── clean_test_results.ps1/.sh # ★ 一键清理脚本 (安全清理历史所有测试目录和临时日志)
│   ├── run_pipeline_all.ps1/.sh # ★ 从头到尾一键主流水线 (MD模拟 -> 统计 -> SVG/PNG/PDF绘图)
│   ├── run_all.ps1 / .sh        # 全流程 MD 主脚本 (支持 100 ns / 测试 / OnlyMD 续跑)
│   ├── run_analysis.ps1 / .sh   # ★ 全流程自动化分析主脚本 (一键执行全部分析与作图)
│   └── analysis/                # 全部结果分析与作图子脚本
├── SCI_Methods_GROMACS_MD_Simulation.docx  # ★ 参照论文对标撰写的完整 SCI 方法学 Word 文档
├── SCI_Methods_GROMACS_MD_Simulation.md    # ★ SCI 方法学对标 Markdown 文档
└── README.md
```

## 2. 运行流程

### 2.1 分子对接（已完成）
论文使用 rhAChE（PDB **4ey6**, 2.40 Å）+ Aβ NMR 结构（PDB **1aml**）在 RosettaDock 上对接。
你已经完成了对接并得到 PAS 位点的最优构象（`*_complex.pdb`）。

### 2.2 体系构建
论文 2.2：封端（蛋白各端 + 断裂处 259/262/492/495）→ 质子化 → TIP3P 水 → 截角八面体盒子 → NaCl 中和。

### 2.3 分子动力学
论文 2.3 的流程对应本套脚本：

| 阶段 | mdp | 说明 |
|------|-----|------|
| 能量最小化 | `1_min.mdp` | 2000 步 EM，重原子 3 kcal/mol·Å² 谐波约束 |
| 升温 (NVT) | `2_heat.mdp` | 约束下 ~1 ns 从低温升至 300 K（NVT，脚本做 50→300 K 梯度）|
| 密度平衡 | `3_equil_npt.mdp` | 约束下 1 ns 恒压（Berendsen）平衡 |
| 预平衡 | `4_equil_npt_free.mdp` | 无约束 1 ns 平衡 |
| 产物动力学 | `5_md.mdp` | 1000 ns NPT，300 K / 1 bar，每 0.2 ns 一帧 → 5000 帧 |

> 测试模式（`TESTING=1`）下，`mdp/test/` 中的 NVT(`2_heat`) / NPT(`3`、`4`) / MD(`5`) 均为 **100 步**。

关键参数与论文的一致性：
- **力场**：默认使用标准 Gromacs 内置的 `amber99sb-ildn` 力场（AMBER 体系下的标准等价方案；如需 `amber14sb` 请通过参数指定：`FORCE_FIELD=amber14sb`）
- **水模型**：TIP3P，截角八面体盒子，周期性边界条件
- **截断**：vdW 与静电均为 **12.0 Å**；长程静电用 **PME**
- **控温/控压**：`v-rescale`（Gromacs 中近似 Langevin）+ `Berendsen`（1 bar）
- **约束**：SHAKE/LINCS 约束与 H 相连的键 → 时间步长 **2 fs**；最小化阶段不约束
- **取样**：每 0.2 ns 保存一帧，共 **5000 帧**

### 2.4 运行命令

```bash
# 测试模式 (快速验证流程跑通, NVT/NPT/MD 各 100 步)
cd gromacs_md/scripts
TESTING=1 ./run_all.sh alllhrc
TESTING=1 ./run_all.sh fllhttr
TESTING=1 ./run_all.sh ylsllqr

# 正式模式 (完整 1000 ns 产物动力学)
./run_all.sh alllhrc
./run_all.sh fllhttr
./run_all.sh ylsllqr
```

- `TESTING=1` 时使用 `mdp/test/`（NVT/NPT/MD 均 100 步），适合先跑通流程。
- 默认（不加变量）使用正式 `mdp/`，产物动力学为完整 1000 ns。
- 每个体系会生成 `md_<前缀>/` 工作目录，产物为 `md.xtc`（5000 帧）+ `md.tpr`。

> 注意：`run_all.sh` 中 `pdb2gmx` 使用 `-ter` 交互式询问末端封端方式。要实现论文的封端，
> 请按提示对每一端选择 ACE（N 端封端）/ NME（C 端封端）。若 AChE 存在内部断裂
> （残基 259/262、492/495），需按片段分别封端，或将 `BREAK_RES` 填为对应残基并拆分链。
> `run_all.sh` 顶部 `CHAIN_ACHE` / `CHAIN_PEP` 需与你的 PDB 链 ID 一致。

## 3. 结果分析（复现论文图1–6 / 表1–2）

分析脚本全部放在 `scripts/analysis/`，进入对应体系工作目录（`cd ../md_<前缀>`）后执行。
先运行 `0_make_index.sh` 建立 AChE / 肽 / 骨架索引（**请修改残基范围**）。

| 脚本 | 复现内容 | 论文章节 |
|------|---------|---------|
| `0_make_index.sh` | 建立索引组 | — |
| `1_rmsd_rmsf.sh` | 复合物/AChE/肽 骨架 RMSD 与 RMSF | 3.1 / 图1 |
| `2_rdf.sh` | 肽围绕 AChE 的 RDF（含四等份可靠性检验）| 3.1 / 图2 |
| `3_sasa.sh` | 复合物 SASA 及收敛性 | 3.2 / 图3 |
| `4_secondary_structure.sh` + `dssp_bins.py` | 肽二级结构倾向（每 50 ns）| 3.2 / 图4 |
| `5_hbond.sh` | AChE–肽 与肽内部氢键 | 3.3 |
| `contacts.py` | 非天然接触（图5A/5B、表1）| 3.3 |
| `bridging_waters.py` | 水介导桥连（图6、表2）| 3.4 |
| `plot_all.py` | ★ **一键生成全套论文出版级矢量/位图 (SVG / PNG / PDF)** 与统计表 | 图 1 - 6 / 汇总 2x3 图 |

### 终极一键全流程：从头模拟到生成出版级 SVG / PNG / PDF 与统计表

你只需要在 PowerShell 中使用一键主流水线脚本 **`run_pipeline_all.ps1`**，即可自动连接执行：动力学模拟 -> 多层次指标分析 -> 批量生成矢量图与数据表：
```powershell
# 1. 完整运行：从头自动构建体系、运行 100 ns 模拟并生成全套论文图表
.\run_pipeline_all.ps1 -System alllhrc

# 2. 已跑通平衡，跳过前 10 步平衡，直接跑 100 ns 产物 MD 并生成全套论文图表
.\run_pipeline_all.ps1 -System alllhrc -OnlyMD

# 3. 仅分析绘图模式：已有 md.xtc 轨迹，直接运行 7 大分析与批量图表导出
.\run_pipeline_all.ps1 -System alllhrc -OnlyAnalysis
```

### 仅分析绘图用法示例

**方法一：Windows PowerShell 原生运行（推荐 Windows 环境使用 `.ps1` 脚本）**
```powershell
# 测试模式 (100步轨迹快速出图)
.\run_analysis.ps1 -System alllhrc -Testing

# 正式生产模式
.\run_analysis.ps1 -System alllhrc
```

**方法二：WSL / Git Bash 运行（`.sh` 脚本已自动支持 WSL 调用 Windows python.exe 时的路径自动转换）**
```bash
# 测试模式(适用 100 步测试轨迹)
TESTING=1 ./run_analysis.sh alllhrc

# 正式模式(适用 1000 ns 完整轨迹)
./run_analysis.sh alllhrc
```
也可单独手动调用绘图脚本生成 SVG / PNG / PDF：
```bash
cd ../md_alllhrc
python3 ../scripts/analysis/plot_all.py --dir . --out ./figures
```
也可以进入具体的工作目录手步执行单个脚本：
```bash
cd ../md_alllhrc
bash ../scripts/analysis/0_make_index.sh
bash ../scripts/analysis/1_rmsd_rmsf.sh
python3 ../scripts/analysis/contacts.py -a 1-537 -p 538-579
python3 ../scripts/analysis/bridging_waters.py -a 1-537 -p 538-579
```

> 残基范围 `-a`（AChE）/`-p`（肽）请按你的体系实际编号调整。

## 4. 与论文结果的对应关系

- 图1A 最佳构象 → 你的 `*_complex.pdb`（PAS 位点对接）
- 复合物稳定、AChE 保持稳定、肽移动最大（图1C/E）→ `1_rmsd_rmsf.sh`
- 肽在 3 Å 内围绕 AChE 的概率高于溶剂（图2）→ `2_rdf.sh`
- SASA 前 300 ns 波动后下降（图3）→ `3_sasa.sh`
- 前 50 ns 以 β-转角为主，之后螺旋占主导至 750 ns（图4）→ `4_secondary_structure.sh`
- Asp1、Phe4 为锚定残基，N 端前四个残基及 His14/Gln15 接触最多（图5A）→ `contacts.py`
- 桥连水数目 > 桥连水分子数，多数短寿命（图6）→ `bridging_waters.py`

## 5. 测试建议

由于 1000 ns 全原子模拟计算量很大（数周量级），建议先用**短测试**验证流程：
- 把 `5_md.mdp` 的 `nsteps` 改为如 `500000`（1 ns）快速跑通；或
- 使用论文提及的 200 ns 级别的体系先做加速模拟验证稳定性。

## 6. 常见问题

- **pdb2gmx 报"chain break / 残基缺失"**：AChE 链在 259/262、492/495 处有断裂，需拆分片段并分别封端，参考 `run_all.sh` 中说明。
- **分析时索引组不存在**：先运行 `0_make_index.sh` 并确认分组名（AChE/Peptide/Backbone）。
- **力场定制或报错**：标准 GROMACS 自带 `amber99sb-ildn`，脚本已设其为默认值。如果你本地安装了第三方 `amber14sb.ff`，可在运行命令前加 `FORCE_FIELD=amber14sb`。
- **末端封端选择（-ter）**：脚本默认开启全自动非交互模式（`INTERACTIVE_TER=0`）。如需手动针对 N/C 端逐链选择封端（ACE/NME/None），请在运行前加 `INTERACTIVE_TER=1`。
- **Windows 下无 bash**：使用 WSL，或安装 Git Bash/Cygwin 后运行。
