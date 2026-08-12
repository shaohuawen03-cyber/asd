# Force field — `amber14sb.ff` (ff14SB)

The paper specifies **ff14SB**. GROMACS does not ship ff14SB, so a community
port is vendored here so the pipeline works out of the box with no manual
download.

* Source: <https://github.com/intbio/gromacs_ff> → `amber14sb_OL15.ff`
* Upstream origin: GROMACS user-contributed force fields /
  <https://fch.upol.cz/ff_ol/gromacs.php>, protein part ported by
  Man Hoang Viet, includes the corrected Na⁺ parameters.
* Reference: Maier et al. *ff14SB: Improving the Accuracy of Protein Side Chain
  and Backbone Parameters from ff99SB*, J. Chem. Theory Comput. 2015, 11, 3696.

The directory bundles OL15 nucleic-acid parameters as well; they are unused by
this project (protein + peptide + water + ions only).

## How the scripts find it

`gmx pdb2gmx` picks up any `*.ff` directory sitting in the **current working
directory** before the ones installed in `$GMXDATA/top`. `scripts/build_system.sh`
therefore copies (or symlinks) `forcefield/amber14sb.ff` into the per-system
working directory and calls `pdb2gmx -ff amber14sb`. Nothing has to be installed
into your GROMACS tree, and nothing needs root.

## Capping groups

`aminoacids.rtp` contains `ACE` and `NME` blocks, which is what the capping step
in `scripts/prepare_structure.py` relies on. Since ACE/NME are named terminating
residues, `pdb2gmx` is run with `-ter` and `None/None` selected for every chain
segment — the AMBER `-ter` machinery must not add charged NH3+/COO- groups on top
of a cap.
