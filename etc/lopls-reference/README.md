# L-OPLS / OPLS-AA reference data

Vendored subset of the GROMACS `oplsaa.ff` force-field directory from

    https://github.com/wesbarnett/lopls
    commit 59ffe9819f53598eb96576c75a6417fef52c1755

fetched 2026-10-28 via `gh api repos/wesbarnett/lopls/contents/...`.

This is James Barnett's GROMACS transcription of the L-OPLS modifications to
OPLS-AA. It is a complete, standard `oplsaa.ff` with the Siu et al. L-OPLS
deltas layered in as `#ifdef LOPLS` / `#ifndef LOPLS` blocks in
`ffnonbonded.itp` and `ffbonded.itp`, so it serves as the single reference
for two separate things used in this library:

- The **`#ifndef LOPLS`** (and unconditional) branches are just the stock,
  unmodified OPLS-AA atom types, charges, and bonded parameters — used as the
  "library charge" source for the PAMA backbone/ester region
  (see `etc/report.md` and the per-monomer provenance notes).
- The **`#ifdef LOPLS`** branches are the actual L-OPLS deltas — used to
  retype the PAMA alkyl side chains beyond the first ester-adjacent CH2
  (see `etc/retype_alkyl_tail_lopls.py`).

Only the files needed for atom types/charges/bonded parameters are kept here
(`atomtypes.atp`, `ffnonbonded.itp`, `ffbonded.itp`, `forcefield.itp`); the
water models, amino-acid building blocks, and other GROMACS force-field
scaffolding from the upstream repo are not needed and were not copied.

## Citation

If parameters sourced from this directory end up in a published PAMA/PAMA-
derived force field, cite:

- S. W. I. Siu, K. Pluhackova, and R. A. Böckmann,
  *J. Chem. Theory Comput.* **8**, 1459–1470 (2012).
  https://doi.org/10.1021/ct200908r
- K. Pluhackova, H. Morhenn, L. Lautner, W. Lohstroh, K. S. Nemkovski,
  T. Unruh, and R. A. Böckmann, *J. Phys. Chem. B* (2015).
  https://doi.org/10.1021/acs.jpcb.5b08569
- James W. Barnett, GROMACS transcription (`wesbarnett/lopls`), per the
  upstream README's request; no explicit license is given, used here as
  reference data (public, unmodified physical-parameter transcription) for a
  force-field library, consistent with how the other third-party parameter
  sets in `polyply/data/` are already incorporated under citation.
