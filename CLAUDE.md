# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Polyply generates input files (topologies/`.itp`) and coordinates (`.gro`) for simulating polymers and
biomacromolecules in molecular dynamics. It builds on `vermouth` (from the Martini force-field ecosystem)
for molecule/force-field data structures and file parsing. Three subcommands make up the CLI, exposed via
the `polyply` script (`bin/polyply`):

- `gen_params` (alias `gen_itp`) — build a topology (`.itp`) for a single molecule from a residue
  sequence/graph and force-field files/library, implemented in `polyply/src/gen_itp.py`
  (`gen_params`/`gen_itp`).
- `gen_coords` — build 3D coordinates for a full system (`.top`), implemented in
  `polyply/src/gen_coords.py` (`gen_coords`).
- `gen_seq` — generate a sequence/graph input file (JSON) for complex (e.g. branched, random
  copolymer) sequences, implemented in `polyply/src/gen_seq.py` (`gen_seq`).

The high-level functions these subcommands call are re-exported from `polyply/__init__.py`.

## Commands

Install (editable, with the vermouth dev dependency it requires):
```bash
pip install git+https://github.com/marrink-lab/vermouth-martinize.git#vermouth
pip install --upgrade -e .
pip install -r requirements-tests.txt
```

Run the full test suite with coverage (as CI does):
```bash
coverage run $(which pytest) -vv
coverage report --omit='*/bin/pytest'
```

Run a single test file or test:
```bash
pytest polyply/tests/test_gen_coords.py
pytest polyply/tests/test_gen_coords.py::test_function_name -vv
```

Lint (CI enforces a minimum score; the wrapper script adds `--fail-under`):
```bash
python run_pylint.py --disable=fixme --fail-under=8.0 polyply
python run_pylint.py --disable=fixme --fail-under=9.5 bin/polyply
```

Build the docs (mkdocs; config at `mkdocs.yml`, content under `docs/`):
```bash
pip install -r requirements-docs.txt
mkdocs serve   # or: mkdocs build
```

## Architecture

**Processor pipeline.** Most transformations on a molecule are implemented as subclasses of
`Processor` (`polyply/src/processor.py`), which defines `run_system(system)` iterating over
`system.molecules` and delegating to a subclass's `run_molecule(meta_molecule)`. Key processors:
`MapToMolecule`, `ApplyLinks`, `ApplyModifications`, `GenerateTemplates`, `Backmap`, `AnnotateLigands`,
`BuildSystem`. `gen_params`/`gen_coords` compose these processors in sequence rather than doing the
transformation work directly themselves.

**Two levels of molecule representation.**
- `MetaMolecule` (`polyply/src/meta_molecule.py`) is a residue-level graph: one node per residue/monomer,
  built from a sequence or graph input (`.txt`, `.json`, FASTA, or a CGsmiles string). This is the
  starting representation for `gen_params`.
- The full atomistic/fine-grained `vermouth.molecule.Molecule` is derived from the `MetaMolecule` by
  `MapToMolecule`, which matches each residue node against `vermouth.Block` definitions from the loaded
  force field, then `ApplyLinks` adds bonded interactions across residue boundaries using
  `vermouth.Link` definitions (matched wherever the residue graph has an edge).

**Topology (`polyply/src/topology.py`)** is the system-level container (parsed from a GROMACS-style
`.top` file via `top_parser.py`) holding multiple molecules, each with its own `MetaMolecule`, plus
system-wide settings (box, molecule counts, distance restraints, etc.) used by `gen_coords`.

**Force-field library** lives under `polyply/data/<forcefield_name>/` (e.g. `martini2`, `martini3`,
`2016H66`, `gromos53A6`, `oplsaaLigParGen`, `parmbsc1`, `ol21`). Each force field directory holds `.ff`
files (vermouth's custom format defining blocks/links per residue) and/or raw `.itp` files. Library
selection happens via `-lib <name>` on the CLI, resolved in `polyply/src/load_library.py`.
`LIBRARY.md` documents which polymers/residues are available in which force field and links to the
corresponding data files — update it when adding new library parameters.

**Coordinate generation (`gen_coords.py`)** builds a system geometrically: `GenerateTemplates` builds a
single-residue 3D template per unique residue type, then `BuildSystem`/`random_walk.py` performs a
multiscale random walk to place residues/molecules at a target density or box size, with support for
cyclic molecules, distance restraints, ligand placement near a host molecule (`annotate_ligands.py`),
and a configurable starting node. `backmap.py` handles mapping/scaling fine-grained coordinates back
from the meta-molecule level.

**Parsers** (`polyply_parser.py`, `top_parser.py`, `ff_parser_sub.py`, `build_file_parser.py`,
`simple_seq_parsers.py`) extend vermouth's `.ff`/`.itp`/`.top` reading with polyply-specific sections
(e.g. sequence macros for `gen_seq`, per-molecule build options for `gen_coords`).

Tests mirror `polyply/src/*.py` one-to-one under `polyply/tests/test_*.py`, with shared fixtures in
`polyply/tests/example_fixtures.py` and input data in `polyply/tests/test_data/`.

## This branch vs. `master`

This branch (`2026-10-28-PAMA-parameters`, descended from `2025-12-01-PAMA-parameters`) is **purely
additive data/tooling work** — it does not touch any file under `polyply/src/`, `polyply/tests/`, or
the CLI. `git diff master...HEAD --stat` shows 35 new files / ~13.8k lines, all under `etc/` and
`polyply/data/oplsaaLigParGen/`. No existing library files are modified; everything here is new.

### Purpose

The branch builds a new **oplsaaLigParGen force field entry for PAMA** — poly(alkyl methacrylate)
copolymers — i.e. an all-atom OPLS-AA parameter set for six methacrylate monomers (EH, LA, MMA, OC, ST,
TD) generated via the [LigParGen server](http://zarbi.chem.yale.edu/ligpargen/), plus the backbone link
(bond/angle/dihedral) definitions needed for polyply to connect them into a polymer chain. This is a
library-content contribution, analogous to the other OPLS/GROMOS entries documented in `LIBRARY.md`,
not a change to polyply's code behavior.

### New library content (`polyply/data/oplsaaLigParGen/`)

- `PAMA_{EH,LA,MMA,OC,ST,TD}.oplsaa.LigParGen.ff` — one block definition per monomer: charge-normalized,
  OPLS-atom-type-renumbered topology derived from the raw LigParGen `.itp` output. The `[ atomtypes ]`
  section is stripped from each (GROMACS/vermouth picks atom types up from the shared nonbonded file
  instead, not from the per-monomer block).
- `PAMA.oplsaa.LigParGen_links.ff` — auto-generated `vermouth.Link` entries for every pairwise
  monomer-monomer (and monomer-cap) backbone connection: bonds/angles/dihedrals across the Calpha–Cbeta
  vinyl attachment point, using standard OPLS-AA CT–CT parameters.
- `CH3n.ff` — a charge-neutral methyl chain-cap residue, needed to terminate a PAMA chain without
  leaving a net charge.
- `README.md` — provenance log (source `.itp`/`.pdb` files and the `cp`/`sed` commands used to produce
  each library entry).

Each monomer occupies a dedicated, non-overlapping block of `opls_NNNN` atom types starting at 1000
(see the table in `etc/README.md`), to avoid clashing with standard OPLS-AA (`< 1000`) types or with
other monomers already in the library.

### Supporting non-library content (`etc/`)

`etc/` holds the maintenance scripts, intermediate/source data, and a parametrization report used to
build the above — these are workflow tooling, not part of the installed `polyply` package:

- `generate_pama_links.py` (+ `.ipynb`) — generates `PAMA.oplsaa.LigParGen_links.ff` from the monomer
  PDB files, given each monomer's Calpha/Cbeta backbone atom names.
- `renormalize-charges.py` — makes a raw LigParGen `.itp` exactly charge-neutral (arbitrary-precision
  arithmetic, residual distributed across the first few atoms).
- `renumber-opls-types.py` — remaps a monomer's `opls_N` types into the next free block, given the
  already-integrated library files.
- `update-ffnonbonded.py` — merges a new monomer's `[ atomtypes ]` section into the shared
  `custom-ffnonbonded.itp`, normalizing LigParGen's column format and sorting by OPLS number.
- `visualize_monomers.py` — RDKit-based 2D SVG rendering of the monomer PDB structures (`*.svg` files in
  `etc/PAMA-monomers/`), for visual QA of connectivity/bond-order inference.
- `PAMA-monomers/` — raw LigParGen source files (`*.pdb`, `itp/*.itp`, charge-normalization logs) with
  provenance notes, one subdirectory per monomer's inputs.
- `terminals/CH3.pdb` — capping-group PDB used when generating terminal links.
- `custom-ffnonbonded.itp` — working copy of the merged nonbonded parameter table (source for the
  deployed library, if/when one is published outside this repo).
- `report.md` — structural/numbering analysis of the six monomers (explains why MMA and TD have
  different LigParGen atom-traversal orders than EH/LA/ST, which matters when locating Calpha/Cbeta).
- `README.md` — the step-by-step runbook for adding a new PAMA monomer to the library (copy source
  files → normalize charges → renumber OPLS types → merge nonbonded params → strip atomtypes section →
  identify backbone atoms → regenerate links → deploy → update docs). **Read this before adding or
  modifying any PAMA monomer**, since the steps must be followed in order for the OPLS type ranges and
  link file to stay consistent.

### Working on this branch

- Changes here should stay confined to `etc/` and `polyply/data/oplsaaLigParGen/`; if a change requires
  touching `polyply/src/`, that is a signal the work has grown beyond "add a library entry" and should
  be scoped/discussed separately.
- Follow `etc/README.md`'s workflow order exactly when adding a monomer — the OPLS renumbering step
  depends on scanning *all* previously-integrated `.ff` files, and the links file is fully regenerated
  (not hand-edited) from the monomer set defined at the top of `generate_pama_links.py`.
- `LIBRARY.md` at the repo root is not yet updated for PAMA/oplsaaLigParGen on this branch; update it
  alongside `polyply/data/oplsaaLigParGen/README.md` when this work is finalized.
