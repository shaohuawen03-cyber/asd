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
| **`v2.6-replot-dssp-percent-rg`** | latest | **Replot pack:** DSSP as **% stacked area + occupancy bars**; fig0 has **no peptide RMSD/RMSF**; **complex Rg** before H-bonds; four single replot scripts + one-click |

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
| `rerun_analysis_four.ps1` | **Yes** |
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
**`v2.6-replot-dssp-percent-rg` (HEAD)** is what you use to **draw** the four systems.

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

# already finished:
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

Do **not** use `-OnlyPlot` the first time on fllhttr: it still needs `gmx gyrate` and complex DSSP.
