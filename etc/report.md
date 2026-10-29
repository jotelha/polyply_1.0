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

---

## 6. Backbone hybridization fix + OPLS-AA/L-OPLS re-parametrization (2026-10-28)

Two corrections were applied to all six monomers, superseding the pure
LigParGen CM1A charges discussed in §§1-5 above.

### 6.1 sp2 -> sp3 backbone fix

LigParGen parametrized each monomer as the free vinyl monomer
(`CH2=C(CH3)-COOR`), so the Calpha-Cbeta bond that becomes the polymer
backbone's single bond after polymerization was still encoded as a vinyl
C=C double bond: `b0=0.1340 nm, kb=459403.2` (OPLS `CM=CM`), sp2 "CM"
nonbonded types (sigma=0.355, epsilon=0.317984), two sp2-planarity improper
dihedrals, and C=C no-rotation RB torsions (`58.576, 0, -58.576, ...`)
around that bond. This was identical across all six monomers.

Fixed by `etc/fix_backbone_hybridization.py`, which replaces these with the
generic sp3 OPLS-AA values already used by `generate_pama_links.py` for the
*inter*-residue backbone link (same source: `oplsaa.ff/ffbonded.itp`):
bond CT-CT (0.15290, 224262.4), angles CT-CT-CT (112.700 deg, 488.273),
CT-CT-HC (110.700 deg, 313.800), HC-CT-HC (107.800 deg, 276.144), and
dihedral CT-CT-CT-HC/HC-CT-CT-HC (RB: 0.6276, 1.8828, 0, -2.5104) for the
(always exactly one H outer atom) dihedrals around that bond; the two
sp2-planarity impropers are dropped. Verified: regenerating a trimer with
`polyply gen_params` shows the inter- and intra-residue backbone bond/angle
parameters are now numerically identical, and no `459403`/`58.576`/`10.460`
signature remains anywhere in the assembled topology.

### 6.2 OPLS-AA library charges for the backbone/ester region

Backbone+ester atoms (Calpha, Cbeta, alpha-methyl C+3H, carbonyl C,
carbonyl O, ester O, first ester-adjacent alkoxy C+its H's) were switched
from LigParGen's 1.14*CM1A charges to standard OPLS-AA library charges,
identified by matching this library's existing *bonded*-parameter values
(which already matched standard OPLS-AA entries, confirming LigParGen's
BOSS backend had assigned them by the same analogy) against named entries
in `etc/lopls-reference/oplsaa.ff/atomtypes.atp`:

| role | opls type | charge |
|---|---|---|
| Calpha (quaternary, all-carbon) | `opls_139` ("AA C", analogy per `opls_465`'s own comment: "for R on C=O, use #135-#139") | 0.000 |
| Cbeta (backbone CH2) | `opls_136` | -0.120 |
| Cbeta H (x2) | `opls_140` | +0.060 |
| alpha-methyl C | `opls_135` | -0.180 |
| alpha-methyl H (x3) | `opls_140` | +0.060 |
| ester carbonyl C | `opls_465` ("AA C: esters") | +0.510 |
| carbonyl O (C=O) | `opls_466` ("AA =O: esters") | -0.430 |
| ester O (C-O-C) | `opls_467` ("AA -OR: ester") | -0.330 |
| first alkoxy C, R'=CH3 (MMA only) | `opls_468` ("methoxy C in esters") | +0.160 |
| its H (x3, MMA only) | `opls_469` ("methoxy Hs in esters") | +0.030 |
| first alkoxy C, R'=CH2-... (EH/LA/OC/ST/TD) | `opls_490` ("C(H2OS) ethyl ester") | +0.190 |
| its H (x2) | `opls_140` (generic HC; no dedicated tabulated value found for this specific H in the vendored reference -- documented analogy) | +0.060 |

For MMA the whole backbone+ester fragment sums to exactly 0.000 with these
values (no residual). For the other five, the generic-HC analogy for the
ethyl-ester CH2's hydrogens leaves a consistent +0.060 e residual per
repeat unit, resolved below.

### 6.3 L-OPLS retyping of the alkyl tail

Tail atoms beyond the first ester-adjacent alkoxy carbon were retyped to
L-OPLS (Siu, Pluhackova & Boeckmann, JCTC 2012) via
`etc/retype_alkyl_tail_lopls.py`:

- plain CH3: charge -0.180 -> -0.222; its H: +0.060 -> +0.074 (LJ unchanged)
- plain CH2: charge -0.120 -> -0.148; its H: +0.060 -> +0.074, and LJ
  epsilon 0.12552 -> 0.11000 (sigma unchanged at 0.250)
- the all-carbon CT-CT-CT-CT dihedral fully inside the tail: RB coefficients
  2.9288, -1.4644, 0.2092, -1.6736 -> 0.518787, -0.230192, 0.896807,
  -1.49134
- branch-point/quaternary carbons within the tail (only EH's 2-ethylhexyl
  branch carbon, atom name `C0E`, auto-detected by heavy/H-neighbor count)
  have no tabulated L-OPLS delta and are left on standard OPLS-AA, same
  reasoning as the backbone quaternary carbon (Calpha).

**Deviation from the original plan sketch:** the ester-junction dihedrals
in these files (`CT-C-OS-CT`, `CT-OS-C-HC`, etc.) use OPLS-AA's `C`/`OS`
bond types, confirmed via exact numeric match against
`etc/lopls-reference/oplsaa.ff/ffbonded.itp`. The lopls repo's L-OPLS ester
deltas are only defined for the `C_2`/`O_2` bond-type variant (a different,
newer OPLS ester atom type not used by these monomers), so there is no
L-OPLS value to substitute at the ester junction -- it stays on standard
OPLS-AA. This is a correction to the plan, not an oversight: substituting
numbers from a mismatched bond-type variant would have been wrong.

### 6.4 Final neutralization and verification

After 6.2+6.3, EH carried a small residual (+0.055 e, from its one
branch-point carbon being left on its original, now-inconsistent LigParGen
charge) and LA/OC/ST/TD each carried exactly +0.060 e (from the ethyl-ester
CH2 analogy noted in §6.2). `renormalize-charges.py` was re-run on all six
monomers (same tool, same convention as the original CM1A normalization in
§1-5) to restore exact neutrality; the per-atom correction is on the order
of 0.001-0.002 e, distributed evenly across each repeat unit.

Verified for all six monomers: `polyply gen_params` assembles a trimer
cleanly (no errors/warnings), with no leftover vinyl-bond signatures
anywhere in the output, exact net-zero charge per trimer, and
`generate_pama_links.py`'s output (`PAMA.oplsaa.LigParGen_links.ff`)
unchanged (the inter-residue link never depended on the fixed intra-residue
parameters). No QM dihedral scans or bulk density/viscosity MD validation
were run in this pass (deferred, per
`/home/jotelha/.claude/plans/golden-beaming-hickey.md`).
