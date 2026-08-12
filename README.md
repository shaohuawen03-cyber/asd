# AChE (4EY6) – peptide PAS complex: GROMACS MD pipeline

A complete, runnable GROMACS workflow that reproduces the protocol of
*"Molecular dynamics simulation of the acetylcholinesterase–β-amyloid peptide
complex"* for **your own docked peptide–4EY6 complexes** (`alllhrc`, `fllhttr`,
`ylsllqr`, … docked at the **PAS**).

Everything here has been executed end to end — build → EM → NVT → NPT →
equilibration → MD → analysis → figures — with GROMACS 2024.4. The bundled
smoke test passes **59/59** checks.

---

## TL;DR — test it now

```bash
source /usr/local/gromacs/bin/GMXRC        # your GROMACS

# 1. self-test on the bundled fixture (~2 min, no input needed)
./tests/run_tests.sh

# 2. your own structures
cp /path/to/*_complex.pdb input/
./scripts/run_all.sh --test -nt 8          # 100-step smoke test, all peptides
./scripts/run_all.sh -nt 16 --gpu          # the real 1 µs production runs
```

Per system, manually:

```bash
./scripts/build_system.sh -f input/alllhrc_complex.pdb -d work/alllhrc
./scripts/make_pas_index.sh -d work/alllhrc      # PAS / 344-361 groups
./scripts/run_md.sh   -d work/alllhrc --test     # 100 steps per stage
./scripts/run_md.sh   -d work/alllhrc            # full production
./scripts/analyze.sh  -d work/alllhrc
python3 scripts/plot_results.py -d work/alllhrc
```

## A run finished — is it any good?

```bash
./scripts/check_run.sh -d work/alllhrc --bundle
```

Applies pass/fail criteria to the things that actually go wrong: stage
completion and trajectory length, temperature/pressure/density/energy drift,
**whether the peptide stayed bound**, receptor stability, LINCS warnings, and
periodic-image distance. It prints `[OK]/[WARN]/[FAIL]` per item with a verdict,
and `--bundle` writes a few-hundred-kB tarball of the numbers and logs —
**no trajectories** — that is easy to share.

Full guidance, including which results *look* wrong but are actually what the
paper reports, is in **[`docs/CHECKING_RESULTS.md`](docs/CHECKING_RESULTS.md)**.

---

## What the test mode does

You asked for NVT, NPT and MD to be **100 steps**. That is exactly what
`mdp/test_*.mdp` contains, and it is what `--test` selects:

| stage | production (`mdp/*.mdp`) | test (`mdp/test_*.mdp`) |
|---|---|---|
| `em`    | 2000 steps | 200 steps |
| `nvt`   | 500 000 steps = **1 ns**, 0→300 K ramp | **100 steps** = 0.2 ps, ramp compressed |
| `npt`   | 500 000 steps = **1 ns** | **100 steps** |
| `equil` | 500 000 steps = **1 ns**, unrestrained | **100 steps** |
| `md`    | 500 000 000 steps = **1000 ns (1 µs)** | **100 steps** |

The test files are otherwise **identical** to the production ones — same force
field, cut-offs, thermostat, barostat, constraints — so a green smoke test means
the production run is wired up correctly. Output frequencies are also scaled
down so the short runs still yield frames to analyse.

---

## The protocol, and how it maps onto GROMACS

| Paper | This repo |
|---|---|
| ff14SB | `forcefield/amber14sb.ff` (vendored port, no install needed) |
| TIP3P water | `-water tip3p` |
| truncated octahedron box | `editconf -bt octahedron` |
| NaCl at physiological conc. + neutral | `genion -neutral -conc 0.15` |
| caps on chain ends **and** breaks (259/262, 492/495) | `prepare_structure.py` detects breaks geometrically and adds ACE/NME |
| 2000 steps EM, 3 kcal/mol·Å² on heavy atoms | `mdp/em.mdp` + `make_posre.py` → 1255.2 kJ/mol/nm² |
| no H-constraints during EM only | `constraints = none` in `em.mdp` |
| 1 ns heating 0→300 K, restrained | `mdp/nvt.mdp`, `annealing = single` |
| 1 ns density equilibration, restrained | `mdp/npt.mdp` |
| 1 ns unrestrained equilibration | `mdp/equil.mdp` |
| 1000 ns production, 300 K / 1 bar | `mdp/md.mdp` |
| Langevin thermostat | `integrator = sd` (stochastic dynamics *is* Langevin) |
| Berendsen barostat | `pcoupl = Berendsen` |
| 12 Å cut-offs, PME | `rvdw = rcoulomb = 1.2`, `coulombtype = PME` |
| SHAKE on X–H, 2 fs | `constraints = h-bonds` (LINCS), `dt = 0.002` |
| frame every 0.2 ns → 5000 frames | `nstxout-compressed = 100000` |

### Three deliberate substitutions

1. **SHAKE → LINCS.** GROMACS' SHAKE does not work with the default domain
   decomposition and is slower. LINCS constrains the same bonds, permits the
   same 2 fs step, and is the standard GROMACS choice.
2. **Langevin thermostat → `integrator = sd`.** In GROMACS, Langevin dynamics is
   an *integrator*, not a `tcoupl` option, so `tcoupl` stays `no`; `tau-t` is
   the inverse friction. This matches the paper's thermostat exactly.
3. **Berendsen barostat is kept for fidelity** but does **not** generate a true
   NPT ensemble. For statistically rigorous production, switch `md.mdp` to
   `pcoupl = C-rescale`. This is called out in a comment in the file.

---

## Repository layout

```
mdp/                    em/nvt/npt/equil/md .mdp   + test_*.mdp (100 steps)
scripts/
  prepare_structure.py  clean, detect chain breaks, add ACE/NME caps
  build_system.sh       pdb2gmx → box → solvate → ions → index   (one command)
  make_posre.py         restraints at 3 kcal/mol·Å² → -DPOSRES_HEAVY
  make_index.py         Receptor / Peptide / Complex / backbone groups
  make_pas_index.sh     PAS, PAS_aromatic, Region344_361, L76/W77, CAS
  run_md.sh             the five stages, resumable, --test for 100 steps
  analyze.sh            all of the paper's analyses in one pass
  contacts.py           Fig 5 + Table 1 (7 Å non-native contacts)
  water_bridges.py      Fig 6 + Table 2 (water-mediated bridges)
  plot_results.py       Figures 1–6 as PNGs, with multi-peptide overlay
  run_all.sh            every structure in input/, then a comparison figure
  check_run.sh          post-run health check + shareable report bundle
forcefield/amber14sb.ff ff14SB port (used automatically)
tests/                  run_tests.sh + a synthetic complex fixture
input/                  ← put your *_complex.pdb here
docs/PROTOCOL.md        stage-by-stage reference and troubleshooting
docs/CHECKING_RESULTS.md  how to judge a finished run, and what to send me
```

---

## Structure preparation

`prepare_structure.py` is the part that usually breaks on docked complexes, so
it is deliberately defensive:

* strips waters, ions, the galantamine ligand, NAG/BMA/MAN glycans and other
  HETATM records;
* maps `MSE→MET`, `HID/HIE/HIP→HIS`, `CYX→CYS`, … and fixes `SE→SD`;
* keeps the first altLoc, strips input hydrogens (pdb2gmx rebuilds them);
* **finds chain breaks geometrically** — any C(i)–N(i+1) distance > 2.5 Å — so
  4EY6's gaps at 259/262 and 492/495 are caught without hard-coding;
* caps every segment with **ACE/NME**, built from ideal internal coordinates
  with a torsion scan that avoids steric clashes;
* reports S–S candidates and writes a JSON summary that downstream scripts read.

Because every terminus is a named cap, `pdb2gmx` is driven with `-ter` and
`None` for each prompt — the AMBER `-ter` machinery must not add charged
NH3⁺/COO⁻ on top of a cap. `build_system.sh` counts the prompts from the JSON
and answers them automatically.

> **Check the segment table it prints.** If your peptide is chain B, you should
> see one segment for it. If the receptor shows far more segments than
> expected, your docked PDB has more missing loops than 4EY6 normally does.

---

## Analyses produced

| Output | Paper |
|---|---|
| `rmsd_{Complex,Receptor,Peptide}_backbone.xvg` | Fig 1 A/C/E |
| `rmsf_{Receptor,Peptide}_backbone.xvg` | Fig 1 B/D/F |
| `rdf_peptide_around_receptor.xvg`, `rdf_quarter1-4.xvg` | Fig 2 A/B |
| `sasa_complex.xvg` | Fig 3 |
| `dssp_count.xvg`, `dssp_peptide.dat` | Fig 4 |
| `contacts_peptide_receptor.dat`, `contacts_peptide_intra.dat`, `contact_partners.dat` | Fig 5, Table 1 |
| `water_bridges.dat` | Fig 6, Table 2 |
| `hbond_peptide_{receptor,intra}.xvg` | §3.3 |
| `mindist_peptide_PAS.xvg` | PAS occupancy |

Trajectories are made whole and fitted on the receptor backbone before
analysis, so the peptide never appears to fly away across a periodic boundary.

---

## Requirements

* **GROMACS ≥ 2021** (tested on 2024.4). `gmx dssp` needs ≥ 2023; older builds
  fall back to `gmx do_dssp`.
* **Python ≥ 3.8** — the analysis scripts use only the standard library.
* **matplotlib** — only for `plot_results.py`.

No AmberTools, no MDAnalysis, no root access.

---

## Cost warning

The production run is **1 µs per peptide**. For a ~76 000-atom system that is
roughly 2–6 weeks per peptide on a single modern GPU, times three peptides.

Sensible options:

* run the three peptides concurrently on separate GPUs;
* shorten `md.mdp` to 100–200 ns (`nsteps = 50000000` … `100000000`) — usually
  enough to see whether the peptide stays at the PAS;
* keep `nstxout-compressed = 100000` so frame counts stay comparable.

Always confirm with `--test` before committing to a long run.

---

## Troubleshooting

**`Residue 'XXX' not found in residue topology database`** — an unhandled
HETATM survived cleaning. Add it to `RESIDUE_ALIASES` in
`prepare_structure.py`, or delete it from the input.

**`Atom H... in residue ... not found in rtp entry`** — run with `-ignh`
(already the default in `build_system.sh`).

**`number of atoms in the .tpr does not match`** — you edited the topology after
building. Re-run `build_system.sh` into a clean directory.

**LINCS warnings / the run explodes in NVT** — usually a clash in the docked
pose. Check the minimisation converged (`em.log`) and inspect
`prepared.pdb` around the peptide.

**`Region344_361` is empty** — your receptor numbering differs from mature AChE.
Use `make_pas_index.sh --offset N`.

See `docs/PROTOCOL.md` for more.
