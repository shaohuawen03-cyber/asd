# Version history & how to re-plot the four v1.0 systems

---

## 0. One-sentence map

| You asked | Answer |
|---|---|
| Which version ran the original 100 ns MD? | **`v1.0-gromacs-native-pipeline`** |
| What do I wait for? | The four `gmx mdrun` jobs. Do not stop them. |
| What do I re-run after they finish? | **Analysis only** (`rerun_analysis_four.ps1` / `run_analysis.ps1`). Never the MD launchers. |
| Which code draws the new figures? | Current analysis (`v2.5`), on the **same v1.0 trajectories**. |

---

## 1. `v1.0-gromacs-native-pipeline` = your original full MD

This is the **first complete protocol** in this project (not the later 12-step split-topology).

**Launchers (do not run again on a finished folder):**

```
run_all_four_100ns_formal.ps1
  → run_100ns_formal.ps1
      → run_all.ps1 / run_all.sh     using mdp/100ns/
      → run_analysis.ps1             (old auto-plots; replace later)
```

Single system:

```powershell
.\run_pipeline_all.ps1 -System alllhrc
.\run_100ns_formal.ps1 -System alllhrc
```

**What v1.0 actually does**

1. One `pdb2gmx` on the whole `input/<sys>_complex.pdb` (`amber99sb-ildn`, TIP3P, `-ignh`)
2. `editconf -c -bt triclinic -d 1.0`
3. `solvate` + `genion -neutral -conc 0.15`
4. EM (`1_min.mdp`)
5. NVT 0→300 K, 1.0 ns, POSRES + simulated annealing (`2_heat.mdp`)
6. NPT restrained 1.0 ns (`3_equil_npt.mdp`, Berendsen, `tau_p=5.0`)
7. NPT free 1.0 ns (`4_equil_npt_free.mdp`)
8. Production 100 ns, 2 fs, 300 K / 1 bar, frame every 20 ps → **5000 frames** (`5_md.mdp`)
9. Output: **`md.xtc` + `md.tpr`** (also `equil_free.gro`, `heat.gro`, `complex.gro`)

**Git tag:** `v1.0-gromacs-native-pipeline` → commit `13cfebe`  
(the first tree on this branch that already contains `run_all.sh` + `mdp/100ns/`)

Check a folder (do not delete):

```powershell
dir F:\0wsh\asd\gromacs_md\md_alllhrc\md.xtc, F:\0wsh\asd\gromacs_md\md_alllhrc\equil_free.gro
```

If those two exist, that system is **v1.0**.

---

## 2. Later tags (analysis only — MD physics unchanged)

| Tag | Commit | Role |
|---|---|---|
| **`v1.0-gromacs-native-pipeline`** | `13cfebe` | **Your production MD.** Keep all `md.xtc`. |
| `v2.3-sci-comprehensive-tables-and-clean-rmsf` | `13cfebe` (same tree) | First auto-plot extras on this branch |
| `v2.4-dssp-perframe-and-peptide-phases` | `496632d` | Peptide-only DSSP (too sparse) — do not use for Fig 4 |
| **`v2.5-complex-dssp-rmsd-diagnosis`** | `fe0a38e` | Complex DSSP engine (keep this DSSP logic) |
| **`v2.6-replot-dssp-percent-rg`** | `6aa707a` / `fce03fa` | Replot pack (had a fig4 legend crash) |
| **`v2.6.1-fix-plot-legend-fig0`** | `5bfbd47` | Fixed fig4 `fontsize` crash; first fig0 with stacked DSSP + Rg (H-bonds were dropped from overview) |
| **`v2.6.2-fig0-hbonds-dssp-peptide`** | latest | **Use this to draw.** fig0 is 2×4 and **H-bonds are back** (after Rg). DSSP = occupancy bars + % lines + residue×time map (the old full-complex stack looked flat). Peptide RMSD/RMSF have their own figure. |

`v2.0-user-custom-pipeline` is the **other** protocol (`run_split_md_workflow.ps1`, `md_0_1.xtc`). It is **not** what you used for the four 100 ns jobs.

---

## 3. After the four v1.0 jobs finish: re-draw figures (no new MD)

Old auto-plots (step after `mdrun` in `run_100ns_formal.ps1`) used the early analysis code.  
New figures = **v1.0 `md.xtc` + current `run_analysis.ps1`** (PBC `md_fit.xtc` + complex DSSP).

```powershell
cd F:\0wsh\asd
git pull origin arena/019ff90e-asd
cd gromacs_md\scripts

# only systems whose mdrun has written md.gro (or md_0_1.gro)
.\rerun_analysis_four.ps1
```

One system:

```powershell
.\run_analysis.ps1 -System fllhttr
```

What this does (does **not** touch `md.xtc`):

1. `0_make_index.sh` — groups AChE / Peptide; if needed rebuild  
   `md_fit.xtc` = `whole → nojump → center -pbc mol -ur compact → fit rot+trans`
2. RMSD / RMSF / RDF / SASA / **complex DSSP** / H-bonds / contacts / waters
3. `plot_all.py` → `md_<sys>/figures/`

`-OnlyPlot` is **not** enough for fllhttr / ylsllqr / ache: they still need step 1–2.  
`alllhrc` you already re-analyzed; the four-system script will just refresh it.

---

## 4. SCI methods text for v1.0 (can go in the paper)

All-atom MD of AChE–peptide complexes was performed with GROMACS using AMBER99SB-ILDN and TIP3P. Each complex PDB was processed with a single `pdb2gmx` (`-ignh`). The solute was centered in a triclinic box with a 1.0 nm solute–wall distance, solvated, and neutralized at 0.15 M NaCl. After energy minimization, the system was heated from 10 K to 300 K over 1.0 ns (NVT, heavy-atom position restraints, v-rescale, GROMACS simulated annealing), then equilibrated for 1.0 ns NPT with restraints and 1.0 ns NPT without restraints (Berendsen barostat, 1 bar, τP = 5.0 ps). Production MD was 100 ns NPT (300 K, 1 bar, 2 fs, LINCS on H-bonds, PME, Verlet 1.2 nm). Coordinates were saved every 20 ps (5000 frames). Trajectories were made whole, jump-corrected, centered, and fitted (`gmx trjconv`) before analysis. Backbone RMSD/RMSF, COM RDF, SASA, DSSP (complex and peptide), hydrogen bonds, contacts, and bridging waters were computed with the v2.5 analysis scripts.

---

## 5. Script safety

| Script | After MD finished? |
|---|---|
| `replot_alllhrc.ps1` / `replot_fllhttr.ps1` / `replot_ylsllqr.ps1` / `replot_ache.ps1` | **Yes** |
| `replot_all_four.ps1` / `rerun_analysis_four.ps1` | **Yes** |
| `run_analysis.ps1 -System X` | **Yes** |
| `run_pipeline_all.ps1 -OnlyAnalysis` | **Yes** |
| `run_fix_dssp_peptide.ps1 -System X` | Yes (DSSP + diagnosis + replot only) |
| `run_all_four_100ns_formal.ps1` / `run_100ns_formal.ps1` / `run_all.sh` | **NO** — delete and restart v1.0 MD |
| `run_split_md_workflow.ps1` / `run_all_four_user_workflow.ps1` | **NO** — different protocol + deletes `md_*` |
| `clean_test_results.ps1` | **NO** — deletes all four `md_*` |

---

## 6. Create / list tags

```powershell
git fetch --tags
git tag -l
git log --oneline --decorate -8
```

`v1.0-gromacs-native-pipeline` points at `13cfebe` (MD protocol).  
**`v2.6.2-fig0-hbonds-dssp-peptide` (HEAD)** is what you use to **draw** the four systems.

---

## 7. How to replot the first two finished systems (and the last two later)

DSSP is still the **v2.5 complex** calculation (`ss_complex_*`, ~530 residues).  
Fig 4 is now a **literature-style percentage figure**: stacked α-helix / β-sheet / turn / bend / coil (0–100%) plus a last-20 ns occupancy bar chart.

Fig 0 (2×3): A RMSD (Complex+AChE) · B AChE RMSF · C RDF · D Complex SASA · E DSSP % · F **Rg**.  
No peptide RMSD/RMSF on the overview. H-bonds stay in `fig_hbonds` (computed after Rg).

```powershell
cd F:\0wsh\asd
git pull origin arena/019ff90e-asd
cd gromacs_md\scripts

# already finished (analysis products exist -> scripts auto OnlyPlot, no 5001-frame redo):
.\replot_alllhrc.ps1
.\replot_fllhttr.ps1

# when the last two mdrun jobs write md.gro:
.\replot_ylsllqr.ps1
.\replot_ache.ps1

# or one click for every finished system:
.\replot_all_four.ps1
```

| Script | What |
|---|---|
| `replot_alllhrc.ps1` | only alllhrc |
| `replot_fllhttr.ps1` | only fllhttr |
| `replot_ylsllqr.ps1` | only ylsllqr |
| `replot_ache.ps1` | only ache monomer |
| `replot_all_four.ps1` | all that have `md.xtc` + `md.gro` |

Switches (all five scripts): `-OnlyPlot` force figures only; `-Full` recompute every analysis.

After the v2.6 fig4 crash you do **not** need `-Full`. Default sees `rmsd_complex_bb.xvg` + `ss_complex_frac.xvg` and only re-runs `plot_all.py`.

---

## 8. v2.6.1 — why the four replot scripts "failed"

The four wrappers themselves were fine (they just call `run_analysis.ps1`).  
Analysis for the finished systems **did finish** (H-bonds / contacts / bridging waters all wrote files).  
`plot_all.py` then died at Figure 4:

```
TypeError: legend() got multiple values for keyword argument 'fontsize'
  safe_legend(ax, loc="upper right", ncol=2, fontsize=8)
  -> ax.legend(frameon=False, fontsize=9, **kwargs)
```

Also in that same v2.6 tree (docs said one thing, code did another):

- fig0-E was still 0–1 **fraction line** plots, not literature stacked %
- fig0-F was still H-bonds, not complex **Rg**

v2.6.1:

1. `safe_legend` uses `setdefault("fontsize", 9)` so callers can pass `fontsize=`
2. fig0-E = stacked helix/sheet/turn/bend/coil **0–100%**; fig0-F = complex/AChE **Rg**
3. H-bonds stay in `fig_hbonds` (analysis order still DSSP → Rg → H-bonds)
4. `replot_*.ps1` share `replot_common.ps1`: skip if `mdrun` still running; OnlyPlot if products exist
5. Prefer v1.0 `md.tpr` / `md.xtc` over `md_0_1.*`
6. Smoke test: `scripts/analysis/test_plot_all_smoke.py` (reproduces the legend crash and checks fig0)

Do **not** stop the four `gmx mdrun` jobs. Do **not** re-run MD launchers.

---

## 9. v2.6.2 — H-bonds back on overview; DSSP figure; peptide RMSD/RMSF standalone

Your fig0 after v2.6.1 was missing H-bonds because v2.6.1 **replaced** panel F (H-bonds) with Rg.  
“加在氢键分析前面” means **insert Rg before H-bonds**, not delete H-bonds.

DSSP numbers were already correct (complex helix ~34%, L~536). The stacked 0–100% panel looks “wrong” because AChE secondary structure is stable (±1%). That is expected, not a parser bug. Fig 4 now shows occupancy bars + a residue×time map, which is what most MD papers use.

Peptide RMSD/RMSF are written to `fig_peptide_rmsd_rmsf` (self-fit RMSD, optional ligand-fit line, per-residue RMSF). Not drawn on fig0. No 3-stage labels painted on the curve.

```powershell
cd F:\0wsh\asd
git pull origin arena/019ff90e-asd
cd gromacs_md\scripts
.\replot_alllhrc.ps1
.\replot_fllhttr.ps1
```

Check: `figures\fig0_summary_all.png` (G = H-bonds), `figures\fig4_secondary_structure.png`, `figures\fig_peptide_rmsd_rmsf.png`.

---

## 10. v2.7 — Unified complex-only figures + apo-vs-complex style compare folders

Why the four synced `fig0_summary_all.png` looked different:

1. `md_alllhrc/figures` were **stale**: generated by a pre-v2.6 plot version
   (6-panel layout, DSSP fraction 0–1 lines, no Rg panel) and `gyrate_*.xvg` was
   never produced for alllhrc (its earlier replot died on the v2.6 legend crash).
2. The other three were v2.6.2 layout but every panel **auto-scaled per system**:
   RMSD 0–0.25/0.30, RMSF 0–0.5/0.6/0.7, RDF peak 404/213/172/44, H-bonds 0–8/0–12.
   Same panels, different axes -> visually "not unified".
3. Remaining differences are **real physics** (peptide binding position/dynamics):
   H-bonds 6.2/3.9/2.9/3.8, peptide RMSD 0.13/0.16/0.23/0.16, RDF peak distance
   1.30/1.22/1.66/2.16 nm. Unified axes make these comparable, not identical.

About "why does the complex fig contain a separate ache result": the "AChE BB" /
"AChE Rg" curves were NOT the separate `ache` system — they were the AChE part of
the same complex trajectory. `md_ache` itself is also a complex (6-residue peptide).
Per request, v2.7 removes all AChE-only curves from the complex figures (fig0/fig1/
fig3/fig_hbonds); AChE-part metrics stay in the tables. Cross-system comparisons now
live in dedicated folders `compare_ache_vs_alllhrc|fllhttr|ylsllqr`.

v2.7 changes:

- `plot_all.py`: `--limits-json` (shared y-limits), complex-only curves,
  fig0 panel B = Complex RMSF (1–537), suptitle shows the system name.
- `plot_compare_systems.py`: 6-panel comparison (RMSD / RMSF / SASA / Rg /
  DSSP last-20-ns grouped bars / AChE–peptide H-bonds) + `compare_summary.csv`
  (mean ± SD per system + delta).
- `unified_replot_and_compare.py` + `run_unified_replot.ps1`: one-click pipeline.
  Computes one shared y-range from all four systems (`unified_limits.json`),
  replots all four fig packs, then writes the three compare folders.
  No mdrun, no MD re-run. If `gyrate_complex.xvg` is missing (alllhrc) it runs
  `6_rg.sh` automatically (needs `md_fit.xtc` + gmx on PATH, as before).
- Smoke test now also asserts fig0 has NO "AChE BB"/"AChE Rg" curves.

```powershell
cd F:\0wsh\asd
git pull origin arena/019ff90e-asd
cd gromacs_md\scripts
.\run_unified_replot.ps1
```

Check: `md_*\figures\fig0_summary_all.png` (identical panels/axes, complex-only)
and `gromacs_md\compare_ache_vs_*\fig_compare.png` + `compare_summary.csv`.
