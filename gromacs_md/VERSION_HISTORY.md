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
the same complex trajectory. **`md_ache` itself is NOT a complex**: `ache` is the
standalone AChE protein control (apo, no peptide, `input/ache.pdb`, resid 1-530),
while `alllhrc/fllhttr/ylsllqr` are the three AChE+peptide complexes
(AChE 1-530 + peptide 531-537). The earlier "6-residue peptide" claim here was
wrong and is corrected in v2.7.2.
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

---

## 11. v2.7.1 — RMSF straight-line artifact removed

Cause: `gmx rmsf` restarts residue numbering per chain, so `rmsf_complex_bb.xvg`
is chain A 4-258 / 265-494 / 498-542 (missing residues 259-264, 495-497 are gaps)
then chain B restarts at 1-7. Plotting that as one polyline connected (542, y)
to (1, y) with a long straight diagonal line across the RMSF panel (fig0 B,
fig1 B and every compare figure — hence "everywhere").

Fix: `plot_rmsf_profile()` splits the profile at every residue-numbering
discontinuity (gap > 1 or restart), keeps real numbering for gap segments,
renumbers the peptide chain to continue after chain A (543-549), and draws
each segment separately — a small visual gap instead of a diagonal line.
`rmsf_segments()` + unit test `test_rmsf_profile_splits_at_chain_numbering_restart`
lock the behavior (4/4 smoke tests pass).

No other figure logic changed. Re-run `.\run_unified_replot.ps1` to refresh
fig0/fig1/compare RMSF panels.

---

## 12. v2.7.2 — ache = standalone apo control: H-bond comparison is complexes-only

**Your four systems (this is how the MD was set up):**

| system | composition | peptide? |
|---|---|---|
| `ache` | standalone AChE monomer (`input/ache.pdb`, resid 1-530) | **no (apo control)** |
| `alllhrc` / `fllhttr` / `ylsllqr` | AChE (1-530) + peptide (531-537) | yes |

So AChE–peptide H-bonds, peptide–AChE RDF, peptide RMSD/RMSF/DSSP **only exist
for the three complexes**. `ache` must never appear in any AChE–peptide H-bond
comparison (no curve, no legend entry, no CSV value).

**What was wrong before v2.7.2**

1. `0_make_index.sh` only knew 537/579-residue systems. For the 530-residue apo
   `ache` it hit the "last 7 residues = peptide" fallback, so residues 524-530 of
   AChE itself were labeled `Peptide` -> bogus `hbond_ache_pep.xvg`, `rdf_pep_ache.xvg`
   etc. were computed and `ache` showed a fake "AChE-Peptide" curve in the
   compare figures.
2. `plot_compare_systems.py` drew panel F (H-bonds) for both systems and its
   legend said "ache AChE-Peptide". Also panels A-D labeled the apo control
   "ache Complex BB/SASA/Rg" — wrong, it is not a complex.
3. `contacts.py` / `bridging_waters.py` / `compute_peptide_ss.py` had the same
   "last 7 residues" fallback and could fabricate peptide data for apo.
4. VERSION_HISTORY section 10 wrongly claimed `md_ache` is "also a complex".

**v2.7.2 changes**

- `plot_all.py`: `is_apo()` (index.ndx `[ Peptide ]` group is authoritative,
  peptide-product files as fallback). Apo figures now say "AChE BB / AChE SASA /
  AChE Rg", fig0 suptitle "AChE MD summary (apo control, no peptide, ACHE, 100 ns)",
  RDF panel shows "N/A — apo control (no peptide)", tables use system label
  "AChE" instead of "Complex". The shared `hbond` y-limit (built from
  AChE-peptide counts, 0-10) is NOT applied to intra-AChE counts (~400):
  `fig_hbonds` of the apo control plots "Intra-AChE" on its own scale.
- `plot_compare_systems.py`: panel F only draws systems that have a peptide;
  legend labels are apo-aware ("ache AChE BB (apo)" vs "alllhrc Complex BB");
  suptitle "ache (AChE apo control) vs alllhrc (AChE-peptide complex)";
  `compare_summary.csv` hbond/RDF/peptide-RMSD cells are empty for apo and
  no delta is computed (metric names are now peptide-agnostic:
  Backbone_RMSD/SASA/Rg_last20ns...).
- `0_make_index.sh`: 530 residues -> monomer mode (`AChE`, `AChE_Backbone`
  only, NO Peptide group) and removes stale peptide-dependent products
  (hbond_ache_pep, rdf_pep_ache*, rmsd_pep_bb, ss_pep_*, inter/intra contacts,
  bridging waters) so old bogus files can never leak into new figures.
- `contacts.py` / `bridging_waters.py` / `compute_peptide_ss.py`: 530-residue
  apo -> empty peptide selection -> skip instead of using the C-terminal
  7 residues as a fake peptide.
- Smoke test: 6 tests, two new ones lock the apo behavior (apo fig0 contains
  no "AChE-Peptide"/"Complex BB"; apo vs complex compare figure draws the
  H-bond curve only for the complex and leaves the apo CSV cells empty).

**How to fix your local copy (analysis + figures only, no MD re-run):**

```powershell
cd F:\0wsh\asd
git fetch origin arena/01a03bd4-asd
git checkout origin/arena/01a03bd4-asd -- `
  gromacs_md/scripts/analysis/plot_all.py `
  gromacs_md/scripts/analysis/plot_compare_systems.py `
  gromacs_md/scripts/analysis/test_plot_all_smoke.py `
  gromacs_md/scripts/analysis/0_make_index.sh `
  gromacs_md/scripts/analysis/contacts.py `
  gromacs_md/scripts/analysis/bridging_waters.py `
  gromacs_md/scripts/analysis/compute_peptide_ss.py `
  gromacs_md/scripts/unified_replot_and_compare.py `
  gromacs_md/VERSION_HISTORY.md

cd gromacs_md\scripts
# 1) rebuild md_ache index in MONOMER mode + remove its stale fake-peptide files:
.\replot_ache.ps1 -Full
# 2) unified axes for the 4 systems + the 3 compare folders:
.\run_unified_replot.ps1
```

Check afterwards:

- `gromacs_md\compare_ache_vs_alllhrc\fig_compare.png` — panel F contains only
  the `alllhrc` AChE-Peptide curve; legend has no `ache` H-bond entry;
  panels A-D legend reads "ache AChE ... (apo)".
- `gromacs_md\compare_ache_vs_*\compare_summary.csv` — hbond/RDF/peptide-RMSD
  cells of `ache` are empty, no delta.
- `gromacs_md\md_ache\figures\fig0_summary_all.png` — suptitle "apo control",
  panel G = "Intra-AChE", RDF panel = "N/A — apo control (no peptide)".
- `gromacs_md\md_ache\hbond_ache_pep.xvg` no longer exists (removed by step 1).

---

## 13. v2.7.3 — Rg auto-fill fixed; index built from topology; two data problems found

**Why alllhrc's Rg panel stayed blank (root cause)**

alllhrc was the only system whose `gyrate_*.xvg` had never been produced. The
v2.7 auto-fill (`ensure_gyrate`) called `bash F:\0wsh\...\6_rg.sh`:

1. `shutil.which("bash")` can pick up WSL's `C:\Windows\System32\bash.exe`,
   which cannot open a `F:\...` Windows path -> the fill silently failed.
2. Even in Git Bash, a failed `gmx gyrate` can leave an **empty/header-only**
   `gyrate_complex.xvg` behind. The old code only checked file *existence*,
   so every later run saw "gyrate present" and the panel stayed empty.
3. `plot_all.py` showed the "Rg missing" notice only when BOTH gyrate files
   were absent; with a stray `gyrate_ache.xvg` the panel rendered **blank**.
4. `replot_common.ps1`'s `$HasRg` also only checked existence.

v2.7.3 fixes (all four):

- `ensure_gyrate`: a file with <2 numeric rows counts as missing and is
  removed; GROMACS is called **directly** (`gmx.exe`/`gmx` on PATH) with
  group fallbacks Protein -> System -> Backbone -> AChE, output validated;
  bash `6_rg.sh` is only a last resort (path converted to `/mnt/<drive>/`
  for WSL bash); failures print the real gmx/bash error tail.
- `6_rg.sh`: same group fallbacks + content validation; stale empty files
  removed before each attempt; works with or without index.ndx.
- `plot_all.py`: Rg panel NEVER silently blank — "Rg missing — run 6_rg.sh"
  notice is shown whenever complex Rg data is absent (even if
  gyrate_ache.xvg exists). Added `xvg_has_data()`.
- `replot_common.ps1`: `$HasRg` requires >=2 data rows; after running
  6_rg.sh it re-verifies and prints a clear warning.

**Index groups now come from the topology, not residue numbers**

The four PDBs are numbered chain A = 4-542 (or 1-542 for ache) + peptide
chain B = 1-7. The old `0_make_index.sh` selections `ri 1-530` /
`ri 531-537` were **wrong for this numbering** (they selected AChE's own
residues 531-537, and cut AChE 531-542 out of the AChE group). v2.7.3
builds the groups from the pdb2gmx itp files instead:

- chain A atoms = 1..NA, chain B atoms = NA+1..NA+NB
  (NA/NB = [ atoms ] counts of topol_Protein_chain_A/B.itp;
   verified against the real data: 8145 / 122-133 atoms).
- No chain B itp (or <2 atoms) -> monomer mode (true apo control).
- `is_apo()` in plot_all now checks the **topology first** (topol.top +
  chain B itp), then index.ndx, then peptide-product files — so a stale
  index can never mislabel a system again.

**Two data problems found in the synced results (2026-08-26, results-sync)**

1. **`md_ache` is NOT the standalone apo control.** Its topology
   (`topol_Protein_chain_B.itp`, built from `complex_clean.pdb`) contains a
   7-residue peptide: ALA-LEU-LEU-LEU-HIS-ARG-CYS — the same ALLLHRC peptide
   as `md_alllhrc` (chain A renumbered 1-542). So all four systems are
   complexes; there is currently **no apo control** in the data. The
   pipeline now decides apo vs complex from the topology, so `ache` will be
   treated as a complex (and appear in H-bond comparisons) until a real apo
   run replaces it (a true apo has no chain B itp and is then excluded
   automatically).
2. **fig5/fig6 data (contacts / bridging waters) were computed on the wrong
   "peptide".** `contacts.py` / `bridging_waters.py` used the residue-range
   guess (531-537), which under the real numbering is AChE's own residues
   531-537. The synced `inter_contacts.csv`, `intra_contacts.csv`,
   `frequent_contacts.tsv`, `bridging_per_residue.csv` are therefore
   AChE-internal, not peptide data. Both scripts now select the peptide by
   chain (`segid B or chainID B`, AChE = protein minus peptide), numbering
   independent. **Re-run contacts + bridging to replace those files**
   (`replot_<system>.ps1 -Full` does it, or run_analysis.ps1).

**What was synced into this branch**

The small result files (xvg / csv / dat / tsv / ndx / log / mdp / top / itp
+ figures for fllhttr / ylsllqr / ache + mdp_templates/) from the user's
results-sync push were merged into `arena/01a03bd4-asd`. `md_alllhrc/figures`
kept the newer branch version. Figures are pre-v2.7.3 output and will be
refreshed by the re-plot below.

**How to fix your local copy**

```powershell
cd F:\0wsh\asd
git fetch origin arena/01a03bd4-asd
git checkout origin/arena/01a03bd4-asd -- `
  gromacs_md/scripts/analysis/plot_all.py `
  gromacs_md/scripts/analysis/plot_compare_systems.py `
  gromacs_md/scripts/analysis/test_plot_all_smoke.py `
  gromacs_md/scripts/analysis/0_make_index.sh `
  gromacs_md/scripts/analysis/6_rg.sh `
  gromacs_md/scripts/analysis/contacts.py `
  gromacs_md/scripts/analysis/bridging_waters.py `
  gromacs_md/scripts/analysis/compute_peptide_ss.py `
  gromacs_md/scripts/unified_replot_and_compare.py `
  gromacs_md/scripts/replot_common.ps1 `
  gromacs_md/VERSION_HISTORY.md

cd gromacs_md\scripts
# 1) rebuild indexes from topology + recompute the analyses whose old
#    outputs were wrong (contacts / bridging) + fill alllhrc Rg:
.\replot_alllhrc.ps1 -Full
.\replot_fllhttr.ps1 -Full
.\replot_ylsllqr.ps1 -Full
.\replot_ache.ps1 -Full
# 2) unified axes + compare folders (auto-fills alllhrc Rg if still missing):
.\run_unified_replot.ps1
```

Check afterwards:

- `gromacs_md\md_alllhrc\gyrate_complex.xvg` exists with ~5000 data rows
  and fig0 panel F shows the Complex Rg curve.
- `md_alllhrc\inter_contacts.csv` now lists peptide residues 1-7
  (chain B), not 531-537.
- Compare figures: `ache` behaves according to its topology — with the
  current data (ache = complex) panel F contains both systems' AChE-peptide
  H-bond curves; once a real apo run replaces md_ache, ache disappears from
  that panel automatically.

**Re-running the true apo control (md_ache)**

New launcher `run_apo_ache_100ns.ps1` does it end-to-end with safety checks:

1. extracts/validates `input/ache.pdb` (single chain, no chain B, ~530
   residues; refuses to run if `input/ache_complex.pdb` exists or the PDB
   has a second chain);
2. renames the current `md_ache` to `md_ache_complex_backup` (instant,
   nothing deleted);
3. runs the same v1.0 100 ns protocol (mdp/100ns) via run_all.ps1;
4. verifies the new topology has no `topol_Protein_chain_B.itp`;
5. re-runs analysis (monomer index) + `run_unified_replot.ps1`, so fig0
   says "apo control" and the compare H-bond panels exclude ache
   automatically.

```powershell
cd F:\0wsh\asd\gromacs_md\scripts
.\run_apo_ache_100ns.ps1            # 正式 100 ns
# .\run_apo_ache_100ns.ps1 -Testing # 先用 5000 步验证贯通
```

After it finishes: `md_ache\figures\fig0_summary_all.png` shows
"AChE MD summary (apo control...)", and
`compare_ache_vs_*\fig_compare.png` panel F contains only the complex's
AChE-Peptide curve. The old complex results stay in
`md_ache_complex_backup\` (analysis scripts ignore it).
