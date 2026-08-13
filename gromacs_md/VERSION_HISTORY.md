# Version history

Two **different MD protocols** live in this repo. They are not the same.

---

## 1. Which protocol is which (how to tell on disk)

| | **A. Original complete pipeline (first version)** | **B. Later 12-step user split-topology** |
|---|---|---|
| **Name in old notes** | `v1.0-gromacs-native-pipeline` | `v2.0-user-custom-pipeline` |
| **Entry scripts** | `run_all.sh` / `run_all.ps1` | `run_split_md_workflow.ps1` |
| **One-click 4 systems** | `run_all_four_100ns_formal.ps1` | `run_all_four_user_workflow.ps1` |
| **Wrapper** | `run_pipeline_all.ps1` → `run_all.ps1` + `run_analysis.ps1` | 12-step then `run_analysis.ps1` |
| **MDP** | `gromacs_md/mdp/100ns/` | `mdp_templates/` or same 100ns if templates missing |
| **How topology is built** | **One** `pdb2gmx` on the whole `*_complex.pdb` | Split `pro.pdb` + `tai.pdb`, two `pdb2gmx`, merge gro |
| **Production files** | **`md.xtc` + `md.tpr`** | **`md_0_1.xtc` + `md_0_1.tpr`** |
| **Other fingerprints** | `equil_free.gro`, `heat.gro`, `complex.gro`, `box.gro` | `pro.pdb`, `tai.pdb`, `pro_processed.gro`, `tai.gro` |

**Look in `md_alllhrc/` (do not delete anything):**

```powershell
dir F:\0wsh\asd\gromacs_md\md_alllhrc\md.xtc, F:\0wsh\asd\gromacs_md\md_alllhrc\md_0_1.xtc, F:\0wsh\asd\gromacs_md\md_alllhrc\equil_free.gro, F:\0wsh\asd\gromacs_md\md_alllhrc\pro.pdb
```

- If you see **`md.xtc` + `equil_free.gro`** → you ran **protocol A** (original complete pipeline). That is the first version.
- If you see **`md_0_1.xtc` + `pro.pdb`/`tai.pdb`** → that folder was (re)built with protocol B.

Analysis scripts already accept **both** `md.xtc` and `md_0_1.xtc`.

---

## 2. Your original 100 ns MD = protocol A

The first complete workflow in this project is:

```
run_all_four_100ns_formal.ps1
  → run_100ns_formal.ps1
      → run_all.ps1 / run_all.sh     (mdp/100ns)
          pdb2gmx (whole complex)
          editconf -bt triclinic -d 1.0
          solvate + genion 0.15 M
          EM → NVT anneal 0→300 K → NPT restrained → NPT free
          mdrun -deffnm md            → md.xtc / md.tpr
      → run_analysis.ps1
```

Single-system equivalent:

```powershell
.\run_pipeline_all.ps1 -System alllhrc
# or
.\run_100ns_formal.ps1 -System alllhrc
```

**Do not re-run those scripts** on a finished system: `run_100ns_formal.ps1` and `run_all.sh` **delete** the work directory / intermediate files and start MD from zero.

The 12-step split-topology scripts were added later as an alternate path. They are **not** the original production protocol.

---

## 3. After the four `mdrun` jobs finish

1. Leave running jobs alone.
2. Pull current **analysis** code (does not change the trajectories).
3. Re-plot only — never call `run_all_four_100ns_formal.ps1` / `run_100ns_formal.ps1` / `clean_test_results.ps1` again.

```powershell
cd F:\0wsh\asd
git pull origin arena/019ff90e-asd
cd gromacs_md\scripts
.\rerun_analysis_four.ps1
```

Or one finished system:

```powershell
.\run_analysis.ps1 -System fllhttr
```

Step 0 of analysis rebuilds `md_fit.xtc` (PBC: whole → nojump → center → fit) from `md.xtc` or `md_0_1.xtc`, then RMSD / complex DSSP / figures with **v2.5** scripts.

`alllhrc` you already re-analyzed; the command will just refresh it.

---

## 4. Git tags (analysis code, not a second MD)

| Tag | Commit | What it is |
|---|---|---|
| `v2.3-sci-comprehensive-tables-and-clean-rmsf` | `13cfebe` | First snapshot on this branch: **protocol A + B both present**; first auto-plots |
| `v2.4-dssp-perframe-and-peptide-phases` | `496632d` | Peptide-only DSSP (too sparse) |
| `v2.5-complex-dssp-rmsd-diagnosis` | latest | Complex DSSP + RMSD diagnosis — **use this to re-plot** |

MD physics (mdp/100ns, `run_all.sh`) did not change between these tags. Only post-processing did.

---

## 5. Script safety

| Script | After MD is done? |
|---|---|
| `rerun_analysis_four.ps1` | **Yes** — analyze finished systems, never deletes |
| `run_analysis.ps1 -System X` | **Yes** |
| `run_fix_dssp_peptide.ps1 -System X` | Yes (DSSP + diagnosis + replot only) |
| `run_pipeline_all.ps1 -OnlyAnalysis` | Yes |
| `run_all.sh` / `run_100ns_formal.ps1` / `run_all_four_100ns_formal.ps1` | **NO** — wipe and restart MD |
| `run_split_md_workflow.ps1` / `run_all_four_user_workflow.ps1` | **NO** — different protocol + deletes `md_*` |
| `clean_test_results.ps1` | **NO** — deletes all four `md_*` |
