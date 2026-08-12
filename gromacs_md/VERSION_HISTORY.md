# GROMACS MD 模拟与多层次轨迹分析项目：版本迭代与 Tag 映射全记录
# Version History & Git Tag Reference Manual

本文档专门用于记录项目开发过程中的所有功能分支版本、Git Tag 标记、代码变更日志以及不同运行流程（**主推系统原生标准化流程** vs **用户自定义 12 步模板流程**）的切换指引，防止在长周期科研与多体系测试中产生版本混淆。

---

## 一、Git Tag 版本里程碑一览表

| Tag 标签名 | 版本定义与功能主题 | 对应 Git Commit 核心说明 | 推荐适用场景 |
|---|---|---|---|
| **`v1.0-gromacs-native-pipeline`** | **主推：GROMACS 标准化显式显卡加速与连续线性退火主流程** | `4e0c02b` / `d75fe6b`<br>• 三斜盒子优化 (`-bt triclinic -d 1.0 nm` 节省 ~35% 溶剂体积)<br>• GROMACS 原生模拟退火 (`annealing=single single`, 1.0 ns 连续线性 0->300 K)<br>• 黄金三步法去 PBC 净化 (`-pbc nojump` -> `-center -ur compact` -> `-fit rot+trans`)<br>• `tau_P = 5.0 ps` 温和恒压平衡消除密度警告 | **★ 推荐作为默认与主要生产流程**<br>适合对 `alllhrc` / `fllhttr` / `ylsllqr` / `ache` 单体等任何体系从零执行高标准 SCI 发表级 100 ns / 1000 ns 模拟。 |
| **`v2.0-user-custom-pipeline`** | **用户保留：12 步自定义 MDP 模板与多体系顺序自动化套件** | `4dc8177`<br>• 严格集成用户个人 12 步工作流：`pdb2gmx` -> `merge_gro` -> `editconf -bt triclinic` -> `ions.mdp` -> `em.mdp` -> `nvt.mdp` -> `npt.mdp` -> `md.mdp (md_0_1)`<br>• 集成用户本地 `mdp_templates/` 模板库与复合物 Top 拓扑合成<br>• 顺次自动化运行：一个跑完出图马上看，脚本后台继续跑下一个体系 | **★ 作为用户自定义定制备用流程**<br>用于重现或兼容用户原有通过个人 `mdp_templates` 执行的动力学模拟方案。 |
| **`v2.1-clean-comparison-plots`** | **去冗余图表精细化版（仅展示复合物整体与受体主干）** | `195dc90`<br>• 优化 2×3 大总图 (`fig0_summary_all`) 绘图规范：图 A 仅重点展示 `Complex BB` 与 `AChE BB` 稳定收敛轨迹；图 B 仅绘制 `AChE BB` 残基柔性分布，消除单体对齐产生的直线与冗余线条<br>• 修复并统一所有测试断点清场安全防线 | 优化 SCI 投稿图表视觉呈现，消除冗余干扰线。 |
| **`v2.2-robust-pbc-and-integrity-auto-pipeline`** | **多源轨迹自适应支持、肽共价检查智能路径定位与 12 步全自动分析闭环** | `79246e3`<br>• **彻底解决 `verify_peptide_integrity.py` 路径找不到 (`md.tpr, md.xtc`) 问题**：支持在 `scripts/` 目录运行自动侦测 `../md_<sys>` 或直接通过 `--sys alllhrc` / `-d ../md_alllhrc` 传入，并在找不到文件时输出清晰的纯中文引导与原因分析<br>• **打通 `run_split_md_workflow` 15 步全自动闭环**：在 12 步跑出 `md_0_1` 后，自动续接第 13 步黄金四步法去 PBC 叠合 (`-pbc whole` -> `-pbc nojump` -> `-center -ur compact` -> `-fit rot+trans`)、第 14 步骨架共价键完整性验证、第 15 步自动生成全体论文出版级图表<br>• **全分析脚本多文件名自适应**：全部 8 个分析与绘图脚本（`0_make_index` ~ `5_hbond`、`contacts.py`、`bridging_waters.py`、`run_analysis`）均无缝兼容 `md_0_1.tpr/xtc` 与 `md.tpr/xtc` | 彻底解决 PowerShell 用户终端留在 `scripts/` 目录检测产物时的路径报错，实现一键从零跑完到出图验证全流程。 |
| **`v2.3-sci-comprehensive-tables-and-clean-rmsf`** | **当前最新：消除 RMSF 假象直线、DSSP 动态容错解析与 SCI 详尽多列指标报告表** | 当前提交<br>• **消除 RMSF 图中跨链对齐直线**：优化 `plot_all.py` 绘图逻辑，仅展示受体主干 `AChE BB` (残基 1-530) 与对齐后的肽主干 `Peptide BB` (残基 531-537)，排除 GROMACS 重名残基合成时产生的 2 点相连假象直线<br>• **DSSP 2025.1 输出解析加强**：放宽 `parse_dssp_num.py` 字段检测，自适应兼容 8~10 列输出形式，确保二级结构数据表正确生成<br>• **SCI 出版级综合统计表全面扩充**：在生成图表后自动导出 **`SCI_Table1_Comprehensive_MD_Metrics.csv`**（含全体和后 20 ns 平衡区均值±标准差、极值与收敛状态）及紧凑多列的 `summary_metrics_wide.csv`，彻底解决指标稀疏和留空问题 | **★ 当前分支默认工作点 (`HEAD`)**<br>专为解答 SCI 论文审稿问题设计，不干扰正在后台运行的模拟任务。 |

---

## 二、如何随时在不同版本与工作流之间自由切换？

### 1. 回溯并查看某一指定 Tag 版本的代码状态
在 PowerShell 或终端中输入以下指令，即可安全切换只读查看历史版本的代码和配置：
```powershell
# 回到主推的 GROMACS 1.0 标准化显卡加速与退火版本：
git checkout v1.0-gromacs-native-pipeline

# 回到用户定制 2.0 版 (12步个人模板流程)：
git checkout v2.0-user-custom-pipeline
```

### 2. 回到最新的开发工作区 (返回默认分支)
```powershell
# 随时切回咱们共享的在线开发分支
git checkout arena/019fbc93-asd
```

---

## 三、两大并存流程核心运行命令手册

### 方案 A（★ 主推默认流程）：`run_pipeline_all.ps1`
采用标准的 GROMACS 2025 连续模拟退火与标准化 `mdp/100ns/` 参数库，自动清理残留并一体化完成模拟与出图：
```powershell
# 1. 对复合物 alllhrc 跑完 100 ns 生产模拟 + 7大分析 + A-F 出版图：
.\run_pipeline_all.ps1 -System alllhrc

# 2. 自动提取并运行 AChE 无结合肽单体对照组 100 ns：
.\run_pipeline_all.ps1 -System ache
```

### 方案 B（用户保留流程）：`run_all_four_user_workflow.ps1`
遵循用户个人 12 步模板脚本与顺序执行队列（一个跑完出图马上看，后台不间断继续跑下一个）：
```powershell
# 顺序跑完四个体系 (alllhrc -> fllhttr -> ylsllqr -> ache) 并在每个完成后即刻出图：
.\run_all_four_user_workflow.ps1
```

---

## 四、关于单独多肽（Peptide BB）分析取舍说明
- **科学建议**：在复合物体系（`Complex BB`）的轨迹分析中，**完全不需要显式分析或绘制单独小肽（Peptide alone）的 RMSD / RMSF**。
- **原因**：因为 Aβ 小肽作为只有 7 个残基的短多肽，在广阔的溶剂环境中其主要科学价值体现在**它是如何抑制或稳定结合在 AChE 酶的外周阴离子位点（PAS）上的**；
- 将分析重心聚焦在 **复合物整体（Complex BB）** 和 **AChE 酶主干（AChE BB）** 的 RMSD/RMSF，以及界面结合氢键、径向分布函数（RDF）、芳香族残基接触表，能够直接有力地支撑论文结论，并保持投稿图表的极致洁净与高专业度。
