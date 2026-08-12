# `input/` — put your docked complexes here

Drop the best-scoring docked poses in this directory, one PDB per peptide:

```
input/alllhrc_complex.pdb
input/fllhttr_complex.pdb
input/ylsllqr_complex.pdb
```

Then:

```bash
./scripts/run_all.sh --test        # 100-step smoke test on every file
./scripts/run_all.sh               # the real 1 µs runs
```

`run_all.sh` names each working directory after the file with `_complex`
stripped, so the three files above become `work/alllhrc`, `work/fllhttr` and
`work/ylsllqr`.

## What the files should contain

* the **AChE receptor** (4EY6) and the **docked peptide** in *separate chains*;
* the peptide positioned at the **PAS**, as produced by your docking run.

Waters, ions, the galantamine ligand and the N-glycans may stay — they are
removed automatically during preparation. Hydrogens may stay too; they are
stripped and rebuilt consistently with ff14SB.

## Chain identification

The peptide is autodetected as the last chain in the file. If that is wrong:

```bash
./scripts/run_all.sh --peptide-chain B
```

Check the segment table printed by the build step. For the three 7-mers above
you should see a single 7-residue segment for the peptide chain.

*(This directory is otherwise empty on purpose — the PDBs are your data and are
not committed to the repository.)*
