# Protocol reference

Stage-by-stage detail of the pipeline, the reasoning behind each GROMACS
setting, and how to change things safely.

---

## 0. Inputs

Put the docked complexes in `input/`:

```
input/alllhrc_complex.pdb
input/fllhttr_complex.pdb
input/ylsllqr_complex.pdb
```

Each is expected to contain the AChE receptor (4EY6) and one docked peptide, in
**separate chains**. The receptor is usually chain A; the peptide is whatever
chain your docking program wrote. If autodetection guesses wrong, pass
`--peptide-chain B`.

The scripts do **not** care whether waters, the galantamine ligand or the
N-glycans are still present — they are removed during preparation, as the paper
specifies ("the galantamine molecule as well as all water molecules were
removed from the binding site of AChE").

---

## 1. `prepare_structure.py` — cleaning and capping

```bash
python3 scripts/prepare_structure.py input/alllhrc_complex.pdb \
        -o work/alllhrc/prepared.pdb --peptide-chain B
```

| Step | Detail |
|---|---|
| HETATM removal | waters, ions, ligands, NAG/BMA/MAN/FUC glycans |
| residue mapping | `MSE→MET` (with `SE→SD`), `HID/HIE/HIP/HSD/HSE/HSP→HIS`, `CYX/CYM→CYS`, `ASH→ASP`, `GLH→GLU`, `LYN→LYS` |
| altLoc | first conformer kept |
| hydrogens | stripped; `pdb2gmx` rebuilds them consistently with ff14SB |
| incomplete residues | dropped if N, CA or C is missing |
| **chain breaks** | detected where C(i)–N(i+1) > `--break-dist` (2.5 Å) or numbering jumps |
| **capping** | ACE before, NME after every segment |

### Why geometric break detection

The paper caps "both ends of the broken parts (residues 259, 262, 492 and
495)". Those numbers are specific to the 4EY6 entry the authors used. Docked
models are frequently renumbered or have slightly different disordered regions,
so hard-coding them silently produces a wrong topology. Measuring the backbone
distance finds the real gaps in *your* file. The printed segment table shows
exactly what was found:

```
chain A seg 0:   GLU4 .. LEU259   ( 256 res)  caps ACE/NME
chain A seg 1: PRO262 .. GLU492   ( 231 res)  caps ACE/NME
chain A seg 2: SER495 .. THR543   (  49 res)  caps ACE/NME
chain B seg 0:   ALA1 .. CYS7     (   7 res)  caps ACE/NME
```

### Cap geometry

Caps are built with NeRF internal coordinates (Amber ideal values: C–N 1.335 Å,
∠CA–N–C 121.7°, …). The ACE placement scans φ in 15° increments and keeps the
orientation with the largest clearance to nearby atoms; NME places its nitrogen
anti to the carbonyl oxygen. This avoids caps landing inside the protein, which
would blow up minimisation.

### Leaving the peptide's real termini charged

The paper caps everything. If you prefer a zwitterionic peptide:

```bash
--free-termini B      # chain B's real N/C termini stay charged
```

Internal breaks are still capped. Note this changes the net charge, which
`genion` compensates for.

---

## 2. `build_system.sh` — solvated, neutral system

```bash
./scripts/build_system.sh -f input/alllhrc_complex.pdb -d work/alllhrc \
    --peptide-chain B --box-dist 1.2 --conc 0.15
```

1. `prepare_structure.py`
2. `pdb2gmx -ff amber14sb -water tip3p -ter -ignh -merge no`
   — answers `None` to every terminus prompt (count read from the JSON report)
3. `make_posre.py` — restraints to 3 kcal/mol·Å², guard → `POSRES_HEAVY`
4. `editconf -bt octahedron -d 1.2 -c`
5. `solvate -cs spc216.gro` (SPC216 is the standard equilibrated box used to
   fill TIP3P systems; the *model* comes from the topology, not this file)
6. `grompp` + `genion -neutral -conc 0.15 -pname NA -nname CL`
7. `make_index.py`

Key options:

| Option | Default | Meaning |
|---|---|---|
| `--box-dist` | 1.2 nm | solute-to-box distance; ≥ 1.2 nm keeps periodic images beyond the 1.2 nm cut-off |
| `--conc` | 0.15 M | physiological NaCl |
| `--box-type` | `octahedron` | as in the paper; `dodecahedron` is cheaper |
| `--free-termini` | – | chains whose real termini stay charged |

The paper's system has ">76,000 atoms". Yours will land in the same ballpark
with `--box-dist 1.2`; a smaller value shrinks it (and the cost) but risks
self-interaction.

### Why the force field is vendored

GROMACS does not ship ff14SB. `forcefield/amber14sb.ff` is a community port,
copied into the working directory so `pdb2gmx -ff amber14sb` finds it without
touching your GROMACS installation. See `forcefield/README.md`.

---

## 3. Position restraints

The paper restrains protein and peptide **heavy atoms** at **3 kcal/mol·Å²**:

```
3 kcal/mol/Å² × 4.184 kJ/kcal × 100 Å²/nm² = 1255.2 kJ/mol/nm²
```

`pdb2gmx` writes `posre_*.itp` at 1000 kJ/mol/nm² under `#ifdef POSRES`.
`make_posre.py` rewrites the value and retags the guard to `POSRES_HEAVY`, which
is the define used by `em.mdp`, `nvt.mdp` and `npt.mdp`. `POSRES_WATER` is left
alone. The atom selection is pdb2gmx's own — exactly the heavy atoms present in
the input coordinates — which matches the paper's definition.

Whenever restraints are on, `grompp` needs `-r <reference>`; `run_md.sh` adds it
automatically by grepping the mdp for `define = ... POSRES`.

---

## 4. `run_md.sh` — the five stages

| Stage | mdp | Length | Restraints | Ensemble |
|---|---|---|---|---|
| `em`    | `em.mdp`    | 2000 steps | yes | – |
| `nvt`   | `nvt.mdp`   | 1 ns, 0→300 K | yes | NVT |
| `npt`   | `npt.mdp`   | 1 ns | yes | NPT |
| `equil` | `equil.mdp` | 1 ns | **no** | NPT |
| `md`    | `md.mdp`    | 1000 ns | no | NPT |

Velocities carry over through checkpoints (`-t nvt.cpt` etc.), so the
trajectory is continuous. `continuation = yes` after NVT prevents constraints
being re-applied to already-constrained coordinates.

```bash
./scripts/run_md.sh -d work/alllhrc --test        # 100 steps per stage
./scripts/run_md.sh -d work/alllhrc --stage md    # one stage
./scripts/run_md.sh -d work/alllhrc -nt 16 --gpu
```

**Resuming.** If `md.cpt` exists, `mdrun -cpi md.cpt` is added automatically.
To extend a finished run:

```bash
gmx convert-tpr -s md.tpr -extend 500000 -o md.tpr   # +500 ns
gmx mdrun -deffnm md -cpi md.cpt
```

### The unrestrained equilibration matters

The paper reports the Phe4–Trp286 π-stacking being lost *during pre-equilibration*.
That happens in `equil`. Skipping it starts production from an artificially
restrained pose.

---

## 5. `analyze.sh`

```bash
./scripts/analyze.sh -d work/alllhrc
./scripts/analyze.sh -d work/alllhrc --skip-bridges     # bridges are the slow part
./scripts/analyze.sh -d work/alllhrc --bridge-dt 1000   # bridges every 1 ns
```

### PBC treatment (step 0)

```
trjconv -pbc mol -center      # complex whole and centred
trjconv -fit rot+trans        # fitted on the receptor backbone
```

Without this, a peptide that crosses a box face produces a nonsense RMSD spike
and contacts vanish. Everything downstream uses the fitted trajectory.

### Contacts (Fig 5, Table 1)

`contacts.py` counts, per peptide residue, receptor atoms within **7 Å**
(the paper's non-native contact criterion) and reports:

* mean contacts per frame, max, occupancy;
* the receptor residues most often touched (>10 times → Table 1);
* intramolecular contacts inside the peptide, with i, i±1 excluded so backbone
  neighbours do not swamp the signal (Fig 5B).

Atom↔residue mapping comes from `gmx editconf -f md.tpr -o topology.pdb`, which
is written in **global** atom order — unlike `gmx dump`, whose indices are
per-moleculetype and do not line up with index files.

### Water bridges (Fig 6, Table 2)

`water_bridges.py` flags a water as bridging peptide residue R in a frame when
it is simultaneously within the cut-off of R and of the receptor. Consecutive
frames form one continuous bridge; a water that leaves and returns counts twice
— which is why bridge counts exceed distinct-water counts, exactly as the paper
notes.

Reported per residue: distinct bridging waters, number of bridges, longest and
mean lifetime, and the fraction lasting a single frame.

Default cut-off is 0.35 nm (standard heavy-atom H-bond distance). For the
paper's literal 3 Å use `--cutoff 0.3`. The 135° angle criterion is not applied
— the geometric criterion alone is the usual GROMACS convention and is far
cheaper; `gmx hbond` output is provided alongside for the strict definition.

This is the slowest analysis (one `gmx select` pass per residue). On a 5000-frame
µs trajectory use `--bridge-dt 1000` (1 ns sampling) first.

---

## 6. Figures

```bash
python3 scripts/plot_results.py -d work/alllhrc
python3 scripts/plot_results.py --compare work/alllhrc work/fllhttr work/ylsllqr \
        -o work/comparison
```

`--compare` overlays the peptides on every figure — the quickest way to see
which one holds the PAS best.

---

## 7. PAS groups

```bash
./scripts/make_pas_index.sh -d work/alllhrc
```

Appends to `index.ndx` (mature human AChE numbering):

| Group | Residues |
|---|---|
| `PAS` | Tyr72, Asp74, Tyr124, Trp286, Tyr341 |
| `PAS_aromatic` | Tyr72, Tyr124, Trp286, Tyr341 |
| `Region344_361` | 344–361 — the paper's main residence region |
| `Contact_L76_W77` | Leu76, Trp77 — the new patch the paper reports |
| `CAS` | Trp86, Glu202, Ser203, His447 |

Each is intersected with `Receptor`, so a peptide residue sharing a number is
never picked up. If a group comes out empty the script says so explicitly —
use `--offset N` for renumbered constructs.

With `PAS` present, `analyze.sh` writes `mindist_peptide_PAS.xvg`, the direct
test of whether your peptide stays at the PAS.

---

## Changing the production length

1 µs × 3 peptides is expensive. To run 200 ns instead:

```
# mdp/md.mdp
nsteps             = 100000000     ; 200 ns
nstxout-compressed = 100000        ; keep 0.2 ns frames -> 1000 frames
```

Leave `nstxout-compressed` alone so frame spacing stays comparable to the paper.

---

## Validation

```bash
./tests/run_tests.sh           # full smoke test, ~2 min
./tests/run_tests.sh --quick   # stop after the build
KEEP=1 ./tests/run_tests.sh    # keep the scratch dir for inspection
```

59 assertions covering: mdp lengths and cut-offs, HETATM removal, capping, break
detection, system size and neutrality, index groups, the 1255.2 kJ/mol/nm²
restraint value, all five MD stages, every analysis output, and figure
rendering.

The fixture `tests/data/test_complex.pdb` is synthetic — a small receptor with
two engineered chain breaks plus a real 7-mer peptide — so the test runs in
minutes and needs no network access.
