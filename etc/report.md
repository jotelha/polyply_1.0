# PAMA Monomer Parametrisation Report

Analysis of the six OPLS-AA/LigParGen force-field entries for
poly(alkyl methacrylate) (PAMA) monomers in
`polyply/data/oplsaaLigParGen/`.  TD is assessed from
`etc/PAMA-monomers/itp/td-normalized.itp` (library file not yet created)
and cross-checked against an independent second LigParGen run
`~/ws/2026-03-09-tdma-monomer/TD_999253/tmp/TD_999253.itp`.

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
| TD | backbone methyl → Calpha → Cbeta → ester C | C01 | C02 |

MMA is the structural outlier: the atom traversal starts from the **ester**
methyl (COO–CH₃) rather than the backbone methyl, placing the ester-group
atoms earlier in the numbering.  This is the expected consequence of methyl
methacrylate having no alkyl tail.

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

The Cbeta–Calpha–EsterC value (120.7°, k = 711.280) was initially flagged
as an OC anomaly in a draft of this report, but verification against all six
files confirmed it is identical everywhere; the draft error arose from
comparing it with the different Methyl–Calpha–EsterC angle (119.7°).

Improper dihedral force constants are uniform across all monomers:
10.460 kJ mol⁻¹ (vinyl sp² centres) and 43.932 kJ mol⁻¹ (carbonyl oxygen
planarity).  All proper dihedrals use the same OPLS-AA Ryckaert–Bellemans
coefficients for each atom-type quartet.

---

## 3. Differences and anomalies

### 3.1 TD — anomalously low ester carbon charge

The partial charge on the ester carbonyl carbon differs significantly between
TD and all other monomers:

| Monomer | Ester C (e) | Carbonyl O (e) | Ester O (e) | COO net (e) |
|---------|-------------|----------------|-------------|-------------|
| EH      | +0.477      | −0.461         | −0.372      | −0.356      |
| LA      | +0.476      | −0.462         | −0.372      | −0.358      |
| MMA     | +0.472      | −0.457         | −0.366      | −0.351      |
| OC      | +0.480      | −0.463         | −0.368      | −0.351      |
| ST      | +0.471      | −0.470         | −0.379      | −0.378      |
| **TD**  | **+0.302**  | **−0.485**     | **−0.404**  | **−0.587**  |

EH, LA, MMA, OC, and ST cluster tightly (ester C +0.471–+0.480; COO net
−0.351 to −0.378).  TD's ester C charge is **36 % lower** (+0.302) and its
COO group carries nearly **60 % more net negative charge** (−0.587 vs.
≈ −0.36).

### 3.2 TD — anomalously low backbone vinyl carbon charges

The Calpha and Cbeta charges in TD are nearly neutral compared to the other
five monomers, correlated with the ester charge anomaly in §3.1:

| Monomer | Calpha (e) | Cbeta (e) |
|---------|------------|-----------|
| EH      | −0.153     | −0.169    |
| LA      | −0.154     | −0.169    |
| MMA     | −0.155     | −0.167    |
| OC      | −0.157     | −0.188    |
| ST      | −0.143     | −0.165    |
| **TD**  | **−0.031** | **−0.051**|

The five established monomers are tightly clustered (Calpha −0.143 to
−0.157; Cbeta −0.165 to −0.188).  TD's backbone carbons are **≈5× less
negative**, consistent with charge density being redistributed onto the ester
oxygens.

The excess negative charge in TD's ester group and the near-neutral backbone
carbons are compensated across the whole molecule (total charge is exactly
zero in all cases), but the redistribution pattern is inconsistent with the
other five monomers given that all share the same methacrylate functionality.

### 3.3 TD charge anomaly confirmed by an independent LigParGen run

The charges from a second, independently submitted LigParGen job
(TD_999253, different job ID from the original TD_88E777) are:

| Atom | TD_88E777 / td-normalized (e) | TD_999253 (e) | Δ (e) |
|------|-------------------------------|---------------|-------|
| C01 (Calpha) | −0.0306 | −0.0308 | −0.0002 |
| C02 (Cbeta)  | −0.0508 | −0.0509 | −0.0001 |
| C03 (ester C)| +0.3022 | +0.3018 | −0.0004 |
| O04 (ester O)| −0.4043 | −0.4041 | +0.0002 |
| O05 (C=O)    | −0.4849 | −0.4849 |  0.0000 |

The differences are entirely within floating-point precision and the
charge-normalisation correction applied to td-normalized.itp.  The anomalous
charge distribution is **fully reproducible** and therefore reflects the
OPLS-AA charge model as applied by LigParGen to this specific molecular
structure, not a numerical artifact of any single submission.

### 3.4 Minor vinyl-H charge scatter

The hydrogens on Cbeta (the =CH₂ group) show small but consistent
variation across monomers:

| Monomer | H on Cbeta (e) |
|---------|---------------|
| EH      | +0.145        |
| LA      | +0.145        |
| MMA     | +0.145        |
| ST      | +0.150        |
| TD      | +0.135        |

ST is marginally higher than the EH/LA/MMA cluster; TD is marginally lower.
These differences (≤ 0.015 e) are small and likely reflect the different
electronic environments of the alkyl tails propagating slightly through the
conjugated vinyl system.  They do not represent a concern.

### 3.5 MMA — expected differences

The first alkyl carbon attached to the ester oxygen is the ester *methyl* in
MMA (C0B = −0.192 e), compared to a CH₂ group in all other monomers
(+0.014 to +0.023 e).  This is chemically expected: the methyl ester dead-end
CH₃ has three hydrogens drawing electron density back onto the carbon, while
the longer-chain CH₂ is flanked by electron-donating methylenes.

---

## 4. Summary

| Finding | Monomers affected | Expected? |
|---------|------------------|-----------|
| Identical bond parameters for all equivalent chemical groups | All 6 | ✓ Yes |
| Identical angles around Calpha (three distinct angles) | All 6 | ✓ Yes |
| Identical improper/proper dihedral coefficients by atom-type quartet | All 6 | ✓ Yes |
| Ester C charge ≈36 % lower; COO net charge ≈60 % more negative | **TD only** | ✗ Unexpected — confirmed by two independent LigParGen runs |
| Backbone vinyl carbons ≈5× less negative | **TD only** | ✗ Correlated with ester charge anomaly |
| Negative first alkyl C; ester methyl rather than alkyl chain | MMA | ✓ Chemically expected |
| Minor vinyl-H charge scatter ±0.015 e | ST (high), TD (low) | ✓ Within normal LigParGen variability |

---

## 5. Recommendations

**TD charge anomaly (§3.1–3.3):**  The anomalous charge distribution is
reproducible across two independent LigParGen submissions and therefore
cannot be resolved by resubmitting.  The most likely cause is that LigParGen
assigned a different OPLS atom type to the TD ester carbonyl carbon
(opls_803, σ = 3.55 Å, ε = 0.293 kJ mol⁻¹) compared to the other monomers,
resulting in a different charge template.  Recommended actions:

1. Compare the OPLS atom type of the TD ester carbonyl carbon (opls_803 in
   the raw itp, opls_1245 after library renumbering) with those of EH
   (opls_1008), LA (opls_1044), MMA (opls_1089), OC (opls_1107), and ST
   (opls_1143).  If the types differ, check whether the assigned type is
   chemically appropriate for a methacrylate ester carbon.
2. If the type assignment is incorrect, manually reassign the ester carbonyl
   carbon to the same OPLS type used by the other five monomers, recompute
   charges via a fresh QM/MM single point, and re-normalise.
3. As a pragmatic interim measure, constrain the TD ester group charges to
   the mean values from the other five monomers (+0.475, −0.463, −0.372 for
   ester C, C=O, and ester O respectively) and distribute the charge
   difference (≈ +0.46 e total) across the backbone and alkyl chain atoms
   proportionally, then re-normalise.
