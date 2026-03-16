# etc/ — PAMA monomer library maintenance scripts

This directory contains the scripts and source data used to build and extend
the custom OPLS-AA/LigParGen force-field entries for poly(alkyl methacrylate)
(PAMA) monomers in `polyply/data/oplsaaLigParGen/`.

## Directory layout

```
etc/
├── generate_pama_links.py   # generates backbone link definitions (.ipynb version also present)
├── generate_pama_links.ipynb
├── renormalize-charges.py   # neutralises total charge of a repeat unit
├── renumber-opls-types.py   # remaps opls_NNN types to avoid library conflicts
├── PAMA.oplsaa.LigParGen_links.ff   # working copy of the generated links file
├── PAMA-monomers/           # LigParGen source PDBs and ITPs
│   ├── *.pdb                # one per monomer (free vinyl form from LigParGen)
│   ├── README.md            # file provenance / dtool UUIDs
│   └── itp/                 # raw and charge-normalised ITPs
│       ├── *.itp
│       ├── *-normalized.itp
│       ├── *.log
│       └── README.md
└── terminals/               # capping-group PDBs used by generate_pama_links.py
    └── CH3.pdb
```

## Monomer library conventions

### Free-monomer form
Each monomer is parametrised as the **free vinyl monomer**
(e.g. CH₂=C(CH₃)–COOR) using the
[LigParGen server](http://zarbi.chem.yale.edu/ligpargen/).
The vinyl C=C double bond serves as the polymer backbone attachment point.
`generate_pama_links.py` generates the inter-monomer bond definitions that
polyply uses to assemble the repeat-unit chain.

### OPLS atom-type numbering
To avoid conflicts with standard OPLS-AA types (< 1000) and between
monomers, each monomer occupies a dedicated consecutive block starting
at 1000:

| Monomer | opls range  |
|---------|-------------|
| EH      | 1000 – 1035 |
| LA      | 1036 – 1083 |
| MMA     | 1084 – 1098 |
| OC      | 1099 – 1134 |
| ST      | 1135 – 1200 |
| TD      | 1201 – 1251 |

### Backbone connection atoms (Calpha / Cbeta)
The polymer backbone link bonds connect **Calpha** of the preceding repeat
unit to **Cbeta** of the following one.  In the free-monomer PDB, Calpha is
the internal quaternary sp² vinyl carbon (=C(CH₃)(COOR)−) and Cbeta is
the terminal sp² vinyl carbon (=CH₂).

| Monomer | Calpha | Cbeta |
|---------|--------|-------|
| EH      | C04    | C05   |
| LA      | C04    | C05   |
| MMA     | C07    | C08   |
| OC      | C03    | C01   |
| ST      | C04    | C05   |
| TD      | C01    | C02   |

---

## Scripts

### `renormalize-charges.py`

Adjusts atomic partial charges in a LigParGen ITP so that the repeat unit
is **exactly charge-neutral** at 8-decimal-place precision.  Uses
arbitrary-precision arithmetic (`decimal.Decimal`) and distributes any
quantisation residual across the first few atoms.

**Usage:**
```bash
python renormalize-charges.py INPUT.itp OUTPUT-normalized.itp
```

**Example** (TD monomer):
```bash
cd PAMA-monomers/itp
python ../../renormalize-charges.py td.itp td-normalized.itp > td.log 2>&1
```

---

### `renumber-opls-types.py`

Remaps every `opls_N` atom-type label in a LigParGen ITP to a new
consecutive block that does not overlap with any already-integrated library
file.  Pass all existing library `.ff` files as trailing arguments; the
script finds the highest `opls_N` in use and assigns new numbers starting
at `max + 1`.  If no existing files are given, numbering starts at 1000
(the library convention).

The replacement is done in a **single regex pass**, so there is no risk of
double-substitution even when old and new number ranges overlap.

**Usage:**
```bash
python renumber-opls-types.py INPUT.itp OUTPUT.ff [EXISTING1.ff ...]
```

**Example** (adding TD after ST):
```bash
cd /path/to/polyply_1.0

python etc/renumber-opls-types.py \
    etc/PAMA-monomers/itp/td-normalized.itp \
    polyply/data/oplsaaLigParGen/PAMA_TD.oplsaa.LigParGen.ff \
    polyply/data/oplsaaLigParGen/PAMA_EH.oplsaa.LigParGen.ff \
    polyply/data/oplsaaLigParGen/PAMA_LA.oplsaa.LigParGen.ff \
    polyply/data/oplsaaLigParGen/PAMA_MMA.oplsaa.LigParGen.ff \
    polyply/data/oplsaaLigParGen/PAMA_OC.oplsaa.LigParGen.ff \
    polyply/data/oplsaaLigParGen/PAMA_ST.oplsaa.LigParGen.ff
```

Then strip the `[ atomtypes ]` section (GROMACS picks these up from the
`.ff` file at simulation time, so they must not appear in the monomer
topology block that polyply reads):

```bash
sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' \
    polyply/data/oplsaaLigParGen/PAMA_TD.oplsaa.LigParGen.ff
```

---

### `generate_pama_links.py` / `generate_pama_links.ipynb`

Generates the `PAMA.oplsaa.LigParGen_links.ff` file that defines all
backbone bond, angle, and dihedral **link** entries for every pair of
monomers (including self-pairs and pairs with the `CH3` / `CH3n` chain
caps).  Inter-monomer parameters default to standard OPLS-AA values for
CT–CT bonds, angles, and torsions.

The script reads the monomer PDB files from `PAMA-monomers/` to identify
connectivity around the backbone attachment atoms.

**Run from this directory (`etc/`):**
```bash
cd etc/
python generate_pama_links.py
# or open and run generate_pama_links.ipynb
```

The output `PAMA.oplsaa.LigParGen_links.ff` is written into `etc/` and
must then be copied to the library:

```bash
cp etc/PAMA.oplsaa.LigParGen_links.ff \
   polyply/data/oplsaaLigParGen/PAMA.oplsaa.LigParGen_links.ff
```

**To add a new monomer**, edit the four variables near the top of the script:

```python
monomer_set = ('EH', 'LA', 'MMA', 'OC', 'ST', 'NEWMON')

monomer_pdb_dict = {
    ...
    'NEWMON': './PAMA-monomers/newmon.pdb',
}

Calpha = {
    ...
    'NEWMON': 'CXX',   # internal quaternary sp2 vinyl carbon
}

Cbeta = {
    ...
    'NEWMON': 'CYY',   # terminal sp2 vinyl CH2
}
```

---

## Full workflow for adding a new monomer

The steps below assume the new monomer is called `NEWMON` and the
LigParGen server has produced `NEWMON.pdb` and `NEWMON.itp`.

### 1. Obtain and store source files

```bash
cp /path/to/NEWMON.pdb  etc/PAMA-monomers/newmon.pdb
cp /path/to/NEWMON.itp  etc/PAMA-monomers/itp/newmon.itp
```

Document provenance in `etc/PAMA-monomers/README.md`.

### 2. Normalise charges

```bash
cd etc/PAMA-monomers/itp
python ../../renormalize-charges.py newmon.itp newmon-normalized.itp \
    > newmon.log 2>&1
```

Document the command in `etc/PAMA-monomers/itp/README.md`.

### 3. Renumber OPLS atom types and create the library entry

Run from the repository root.  Pass **all** existing library files so the
script can determine the next free block:

```bash
python etc/renumber-opls-types.py \
    etc/PAMA-monomers/itp/newmon-normalized.itp \
    polyply/data/oplsaaLigParGen/PAMA_NEWMON.oplsaa.LigParGen.ff \
    polyply/data/oplsaaLigParGen/PAMA_EH.oplsaa.LigParGen.ff \
    polyply/data/oplsaaLigParGen/PAMA_LA.oplsaa.LigParGen.ff \
    polyply/data/oplsaaLigParGen/PAMA_MMA.oplsaa.LigParGen.ff \
    polyply/data/oplsaaLigParGen/PAMA_OC.oplsaa.LigParGen.ff \
    polyply/data/oplsaaLigParGen/PAMA_ST.oplsaa.LigParGen.ff \
    polyply/data/oplsaaLigParGen/PAMA_TD.oplsaa.LigParGen.ff
    # ... add any further already-integrated files here
```

Strip the `[ atomtypes ]` section:

```bash
sed -i '/^\[ atomtypes \]/,/^\[ /{ /^\[ /!d; /^\[ atomtypes \]/d }' \
    polyply/data/oplsaaLigParGen/PAMA_NEWMON.oplsaa.LigParGen.ff
```

Update the OPLS range table in `polyply/data/oplsaaLigParGen/README.md`
and the backbone-atom table above.

### 4. Identify backbone atoms (Calpha / Cbeta)

Inspect the monomer PDB connectivity:

- **Calpha** — the internal sp² vinyl carbon; has bonds to the methyl group,
  the ester carbonyl, and the double bond to Cbeta.  Bond length to Cbeta
  ≈ 0.134 nm.
- **Cbeta** — the terminal sp² vinyl CH₂; has two hydrogens and the double
  bond to Calpha.  Bond length to Calpha ≈ 0.134 nm.

Add both atom names to the `Calpha` and `Cbeta` dicts in
`generate_pama_links.py`.

### 5. Regenerate backbone links

```bash
cd etc/
python generate_pama_links.py
```

### 6. Deploy the updated links file

```bash
cp etc/PAMA.oplsaa.LigParGen_links.ff \
   polyply/data/oplsaaLigParGen/PAMA.oplsaa.LigParGen_links.ff
```

### 7. Update library documentation

- `polyply/data/oplsaaLigParGen/README.md` — add `cp` / `sed` provenance
  commands for the new monomer.
- `polyply/data/oplsaaLigParGen/citations.bib` — verify LigParGen citation
  covers the new `opls_` type range.
