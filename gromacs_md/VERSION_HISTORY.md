# Version history (this repository)

This file records **actual git commits and tags in this clone**.
Earlier chat notes mentioned `v1.0`–`v2.2` (`4e0c02b`, `4dc8177`, …). Those
commits/tags are **not in this repository** (`git log` only has the commits
below). Do not `git checkout v1.0-gromacs-native-pipeline` here — the tag
did not exist until it is created on a real commit.

Create / list tags:

```powershell
git tag -l
git log --oneline --decorate
```

---

## Commits that exist here

| Tag | Commit | Date | What changed for MD / analysis |
|---|---|---|---|
| **`v2.3-sci-comprehensive-tables-and-clean-rmsf`** | `13cfebe` | 2026-08-12 | Initial tree on this branch. RMSF no longer draws the Complex BB diagonal. DSSP parser only relaxed column count (still 50 ns bins → 2 points). SCI Table 1 export. Full GROMACS workflow scripts + `md_alllhrc` figures from a 100 ns run. |
| **`v2.4-dssp-perframe-and-peptide-phases`** | `496632d` | 2026-08-13 | Per-frame peptide DSSP + 1 ns bins + GROMACS 2025 `~`/`=` parse. **Mistake:** still DSSP on Peptide only (6–7 residues, looks empty); time axis used dt=0.1 ns so 5001 frames became fake 500 ns; painted RMSD stairs as a “feature”. |
| **`v2.5-complex-dssp-rmsd-diagnosis`** | `HEAD` (this commit) | 2026-08-13 | **DSSP primary target = complex (Protein / AChE+peptide, ~530 residues).** Peptide DSSP is a supplement. Time axis: 5000 frames → 0–100 ns. RMSD stairs diagnosed (not PBC). PowerShell messages in English (no GBK mojibake). This file lists real tags. |

---

## How to tag these versions (run once after pull)

```powershell
git tag -a v2.3-sci-comprehensive-tables-and-clean-rmsf 13cfebe -m "v2.3 SCI tables and clean RMSF"
git tag -a v2.4-dssp-perframe-and-peptide-phases 496632d -m "v2.4 per-frame peptide DSSP (peptide-only, time-axis bug)"
git tag -a v2.5-complex-dssp-rmsd-diagnosis HEAD -m "v2.5 complex DSSP + RMSD jump diagnosis"
git push origin --tags
```

---

## MD analysis commands for the current tag (`v2.5`)

Only the broken parts (complex DSSP + RMSD diagnosis + replot):

```powershell
cd gromacs_md\scripts
.\run_fix_dssp_peptide.ps1 -System alllhrc
```

Expect:

- `md_alllhrc\ss_complex_summary.txt` — helix/sheet of **AChE**, string length ~500+
- `md_alllhrc\ss_pep_summary.txt` — 7-mer supplement (coil/turn is normal)
- time range **0–100 ns**, not 0–500 ns
- `peptide_rmsd_jump_diagnosis.txt` — why 23/56 ns stairs are not PBC

---

## Peptide vs AChE (they are separable)

| Object | Residues | Index group | RMSD fit group |
|---|---|---|---|
| AChE | 1–530 | `AChE` / `AChE_Backbone` | self |
| Peptide | 531–537 | `Peptide` / `Peptide_Backbone` | self (internal) or fit AChE (pose) |
| Complex | 1–537 | `Backbone` / `Protein` | both together (PBC-sensitive) |

Topology uses split `pro.pdb` + `tai.pdb`. Analysis never treats them as one chain except when you ask for Complex Backbone RMSD.

---

## Workflow scripts (unchanged roles)

- Default MD: `run_pipeline_all.ps1 -System alllhrc`
- User 12-step: `run_all_four_user_workflow.ps1`
- Analysis only: `run_analysis.ps1 -System alllhrc`
- Fix DSSP/RMSD only: `run_fix_dssp_peptide.ps1 -System alllhrc`
