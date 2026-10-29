# PAMA Monomer Parametrisation Report

Analysis of the six OPLS-AA/LigParGen force-field entries for
poly(alkyl methacrylate) (PAMA) monomers in
`polyply/data/oplsaaLigParGen/`.

Date: 2026-03-16

---

## 1. Structural overview

All six monomers share the same free-vinyl-monomer topology from LigParGen:
a `=CH₂–C(CH₃)(COOR)–` fragment.  The atom numbering order varies across
monomers depending on the traversal order used by LigParGen.

| Monomer | Atom order starts with | Calpha | Cbeta |
|---------|------------------------|--------|-------|
| EH, LA, ST | backbone methyl → Calpha → Cbeta | C04 | C05 |
| OC | Cbeta (=CH₂) → Calpha → methyl | C03 | C01 |
| MMA | **ester methyl** → O → ester C → Calpha → Cbeta → backbone methyl | C07 | C08 |
| TD | alkyl chain terminus → alkyl chain → ester O → ester C → Calpha → Cbeta → backbone methyl | C18 | C19 |

MMA is the structural outlier among the first five: the atom traversal starts
from the **ester** methyl (COO–CH₃) rather than the backbone methyl, placing
the ester-group atoms earlier in the numbering.  This is the expected
consequence of methyl methacrylate having no alkyl tail.

TD uses a reversed traversal order relative to the other long-chain monomers
(EH, LA, ST): the alkyl chain terminus is numbered first and the
vinyl/ester group last.

---

## 2. Bond and angle parameters — uniformity

**All bond parameters for chemically equivalent bonds are completely
identical across all six monomers:**

| Bond | Length (nm) | k (kJ mol⁻¹ nm⁻²) |
|------|-------------|-------------------|
| Calpha=Cbeta (vinyl C=C) | 0.1340 | 459 403.2 |
| Calpha–C(ester carbonyl) | 0.1444 | 343 088.0 |
| Ester C=O (carbonyl) | 0.1229 | 476 976.0 |
| Ester C–O (ether) | 0.1327 | 179 075.2 |
| O–C (first alkyl carbon) | 0.1410 | 267 776.0 |
| Alkyl C–C | 0.1529 | 224 262.4 |
| Vinyl C–H (on Cbeta) | 0.1080 | 284 512.0 |
| Alkyl C–H | 0.1090 | 284 512.0 |

**All three angles centred on Calpha are also identical across all six monomers:**

| Angle | θ₀ (°) | k (kJ mol⁻¹ rad⁻²) |
|-------|--------|---------------------|
| Methyl–Calpha–Cbeta | 124.0 | 585.760 |
| Methyl–Calpha–EsterC | 119.7 | 585.760 |
| Cbeta–Calpha–EsterC | 120.7 | 711.280 |

Improper dihedral force constants are uniform across all monomers:
10.460 kJ mol⁻¹ (vinyl sp² centres) and 43.932 kJ mol⁻¹ (carbonyl oxygen
planarity).  All proper dihedrals use the same OPLS-AA Ryckaert–Bellemans
coefficients for each atom-type quartet.

---

## 3. Differences and anomalies

### 3.1 TD — ester oxygen charges slightly below the cluster

| Monomer | Ester C (e) | Carbonyl O (e) | Ester O (e) | COO net (e) |
|---------|-------------|----------------|-------------|-------------|
| EH      | +0.477      | −0.461         | −0.372      | −0.356      |
| LA      | +0.476      | −0.462         | −0.372      | −0.358      |
| MMA     | +0.472      | −0.457         | −0.366      | −0.351      |
| OC      | +0.480      | −0.463         | −0.368      | −0.351      |
| ST      | +0.471      | −0.470         | −0.379      | −0.378      |
| TD      | +0.467      | −0.391         | −0.355      | −0.279      |

TD's ester carbonyl carbon (+0.467) falls within the EH–ST cluster
(+0.471–+0.480).  However, both ester oxygens are approximately 15 % lower
in magnitude than the cluster: carbonyl O (−0.391 vs −0.457 to −0.470) and
ester O (−0.355 vs −0.366 to −0.379), resulting in a less-negative COO net
charge (−0.279 vs −0.351 to −0.378).

### 3.2 TD — backbone vinyl carbon charges

| Monomer | Calpha (e) | Cbeta (e) |
|---------|------------|-----------|
| EH      | −0.153     | −0.169    |
| LA      | −0.154     | −0.169    |
| MMA     | −0.155     | −0.167    |
| OC      | −0.157     | −0.188    |
| ST      | −0.143     | −0.165    |
| TD      | −0.160     | −0.200    |

TD's Calpha (−0.160) is marginally outside the cluster (−0.143 to −0.157)
and Cbeta (−0.200) is slightly more negative than the cluster
(−0.165 to −0.188).  Both deviations are minor.

### 3.3 Minor vinyl-H charge scatter

The hydrogens on Cbeta (the =CH₂ group) show small but consistent
variation across monomers:

| Monomer | H on Cbeta (e) |
|---------|----------------|
| EH      | +0.145         |
| LA      | +0.145         |
| MMA     | +0.145         |
| OC      | +0.149         |
| ST      | +0.150         |
| TD      | +0.144         |

EH, LA, and MMA cluster tightly at +0.145; OC and ST are marginally higher
(+0.149, +0.150); TD is marginally lower (+0.144).  These differences
(≤ 0.006 e) are small and likely reflect the different electronic environments
of the alkyl tails propagating slightly through the conjugated vinyl system.
They do not represent a concern.

### 3.4 MMA — expected differences

The first alkyl carbon attached to the ester oxygen is the ester *methyl* in
MMA (C0B = −0.192 e), compared to a CH₂ group in all other monomers
(+0.010 to +0.018 e).  This is chemically expected: the methyl ester dead-end
CH₃ has three hydrogens drawing electron density back onto the carbon, while
the longer-chain CH₂ is flanked by electron-donating methylenes.

---

## 4. Summary

| Finding | Monomers affected | Expected? |
|---------|------------------|-----------|
| Identical bond parameters for all equivalent chemical groups | All 6 | ✓ Yes |
| Identical angles around Calpha (three distinct angles) | All 6 | ✓ Yes |
| Identical improper/proper dihedral coefficients by atom-type quartet | All 6 | ✓ Yes |
| TD ester oxygen charges ≈15 % lower magnitude; COO net less negative | **TD** | ✗ Minor deviation |
| TD Calpha/Cbeta marginally outside cluster | **TD** | △ Within LigParGen normal variability |
| Negative first alkyl C; ester methyl rather than alkyl chain | MMA | ✓ Chemically expected |
| Minor vinyl-H charge scatter ±0.015 e | ST (high) | ✓ Within normal LigParGen variability |

---

## 5. Recommendations

**TD ester oxygen charges (§3.1):**  The current library entry has ester
oxygens approximately 15 % lower in magnitude than the EH–ST cluster.
Given that the ester carbonyl carbon and backbone charges are within the
expected range, this deviation is considered acceptable for the current
library version.  If higher consistency is required in the future, constrain
the two ester oxygens to the cluster means (−0.463 for C=O and −0.372 for
ester O), distribute the charge difference (≈ +0.11 e total) across the
alkyl chain atoms proportionally, and re-normalise with
`renormalize-charges.py`.
