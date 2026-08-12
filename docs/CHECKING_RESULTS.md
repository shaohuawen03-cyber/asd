# 结果检查与文件同步

第一个体系跑完之后，怎么判断结果对不对，以及要把哪些文件同步给我。

---

## 一、最快的路径

```bash
source /usr/local/gromacs/bin/GMXRC
cd /path/to/repo

# 如果还没做分析，先补上
./scripts/make_pas_index.sh -d work/alllhrc     # PAS / 344-361 分组
./scripts/analyze.sh        -d work/alllhrc

# 体检 + 打包
./scripts/check_run.sh -d work/alllhrc --bundle
```

`check_run.sh` 会逐项打印 `[ OK ] / [WARN] / [FAIL]`，末尾给出结论；
`--bundle` 额外生成一个 **几百 kB** 的 `alllhrc_report.tar.gz`（**不含轨迹**）。

把那个 tar.gz 发给我就够了 —— 见下面第三节。

---

## 二、七项检查分别在看什么

### 1. 各阶段是否跑完 / 产物轨迹有多长

检查 `em/nvt/npt/equil/md` 五个 `.gro` 是否都在，并从 `md.log` 读出性能，
估算总时长。同时数轨迹帧数：

* **> 4000 帧** — 接近论文的 5000 帧，统计可用
* **400–4000 帧** — 可以看趋势，但统计强度不如论文
* **< 400 帧** — 基本是测试运行，**不能用来下结论**

### 2. 热力学量

| 量 | 合格范围 | 不合格意味着 |
|---|---|---|
| 温度 | 297–303 K | 控温没生效 |
| 压力（**平均值**） | −50 ~ +50 bar | 控压异常 |
| 密度 | 970–1050 kg/m³ | 溶剂化或加离子出错 |
| 势能漂移 | < 1 % | 体系在缓慢崩坏 |

> 瞬时压力在几百 bar 范围剧烈涨落是**正常的**，小体系尤其如此。
> 只有**平均值**有意义，别被单帧数值吓到。

### 3. 肽有没有留在蛋白上 ← **最关键**

看 `mindist_peptide_receptor.xvg` 的平均最小距离：

* **< 0.5 nm** — 全程接触，复合物稳定 ✅
* 0.5–1.0 nm — 松散结合
* **> 1.0 nm** — 肽跑了，复合物解离 ❌

再看 `mindist_peptide_PAS.xvg`：

* < 0.5 nm — 占据 PAS
* 0.5–1.2 nm — **在 PAS 附近但不在 PAS 上**

第二种情况**不是 bug**。论文里 Aβ 恰恰就是这样：预平衡阶段 Phe4–Trp286 的
π-π 堆积就丢了，肽随后迁移到**邻近的 344–361 区域**，只是间歇性地回访 PAS
（702–963 ns）。所以肽离开 PAS 但仍贴在表面 344–361 附近，**正是在复现论文结论**。

### 4. 结构稳定性

* **受体骨架 RMSD**：论文报告 AChE < 2 Å。你的应当 < 3 Å；> 5 Å 说明蛋白解体
* **肽骨架 RMSD / RMSF**：论文里 Aβ 的 RMSF 高达 4–14 Å。
  **小肽剧烈涨落是正常的、也是预期的**，不是错误

### 5. 日志里的警告

统计 `LINCS WARNING`、SETTLE 失败、`Fatal error`。
零警告最好；偶发几条 LINCS 警告通常无害，成百上千条说明体系不稳。

### 6. 周期性镜像距离

复合物与自身周期镜像的最小距离应 **> 2.4 nm**（两倍截断）。
低于 1.2 nm 说明盒子太小，静电上蛋白在和自己作用，结果不可信。

### 7. 分析产物

确认 `analysis/` 下的文件都在，并打印接触数最高的几个肽残基
（对应论文 Fig 5A / Table 1）。

---

## 三、要同步给我的文件

### 推荐：用打包

```bash
./scripts/check_run.sh -d work/alllhrc --bundle
```

生成 `work/alllhrc_report.tar.gz`（几百 kB），里面是：

```
check_report.txt        体检结论
prepared.json           清理/封端/断链检测记录
analysis/*.xvg *.dat    全部分析数据
logs/em,nvt,npt,equil   完整日志
logs/md.log.trimmed     md 日志的首尾各 250 行
mdout/*.mdp             grompp 实际使用的完整参数
molecules.txt           拓扑的分子组成
index_groups.txt        索引分组清单
frame_first.pdb         第一帧复合物坐标
frame_last.pdb          最后一帧复合物坐标
```

**同步方式**（任选）：

1. **推到仓库分支**（最省事）
   ```bash
   mkdir -p reports && cp work/alllhrc_report.tar.gz reports/
   git add -f reports/alllhrc_report.tar.gz
   git commit -m "alllhrc run report"
   git push origin arena/019fbc96-asd
   ```
   然后告诉我一声，我会拉下来看。

2. **直接把 `check_report.txt` 的内容贴给我** —— 如果只想让我快速判断，
   这一个文件通常就够。

### 不要发的东西

`md.xtc`、`md.trr`、`*.cpt`、`*.tpr` —— 动辄几个 GB，仓库放不下，
而且我需要的所有数字都已经在分析产物里了。

### 如果只想手动挑几个文件

最小集合（按重要性排序）：

```
work/alllhrc/check_report.txt
work/alllhrc/analysis/mindist_peptide_receptor.xvg
work/alllhrc/analysis/mindist_peptide_PAS.xvg
work/alllhrc/analysis/rmsd_Receptor_backbone.xvg
work/alllhrc/analysis/rmsd_Peptide_backbone.xvg
work/alllhrc/analysis/contacts_peptide_receptor.dat
work/alllhrc/analysis/contact_partners.dat
work/alllhrc/md.log          (可只取首尾各 200 行)
work/alllhrc/prepared.json
```

---

## 四、自己先看图

```bash
python3 scripts/plot_results.py -d work/alllhrc
```

在 `work/alllhrc/analysis/` 下生成 6 张 PNG，对应论文 Fig 1–6。

三个肽都跑完后，横向比较：

```bash
python3 scripts/plot_results.py --compare \
    work/alllhrc work/fllhttr work/ylsllqr -o work/comparison
```

**看图时重点关注：**

* **Fig 1 RMSD** — 受体曲线应在 1–2 Å 处走平；肽的曲线高、有跳变都正常，
  跳变对应结合模式切换（论文在 ~100 ns 和 ~400 ns 各有一次）
* **Fig 2 RDF** — 论文的判据是 3 Å 处 g(r) > 1，且四等分曲线形状一致。
  四条曲线重合 = 采样收敛
* **Fig 3 SASA** — 论文后期呈下降趋势，对应复合物收紧
* **Fig 5 接触数** — 哪些残基在锚定蛋白，对应 Table 1

---

## 五、常见"看起来不对但其实正常"的情况

| 现象 | 判断 |
|---|---|
| 肽离开了 PAS | **正常**。论文中 Aβ 也离开了，主要驻留在 344–361 |
| 肽 RMSD/RMSF 很大（> 5 Å） | **正常**。短肽本来就柔性，论文报告 4–14 Å |
| 瞬时压力 ±500 bar | **正常**。只看平均值 |
| 预平衡后 π-π 堆积消失 | **正常**。论文明确写了这一点 |
| 受体 RMSD > 5 Å | **不正常**。检查封端和断链检测是否正确 |
| 密度 ≠ ~1000 kg/m³ | **不正常**。溶剂化/加离子有问题 |
| 肽最小距离 > 1 nm 且不回来 | **不正常**（或该肽确实不结合，也是有意义的结论） |

---

## 六、如果 `check_run.sh` 报 FAIL

把 `check_report.txt` 发我，同时附上：

* 受体 RMSD 异常 → `prepared.json` + 构建日志里的分段表
* LINCS 爆炸 → `em.log` 和 `nvt.log`
* 密度不对 → `molecules.txt`
* 起步就崩 → `mdout/mdout_em.mdp`
