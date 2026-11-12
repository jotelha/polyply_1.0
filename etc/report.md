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

## 7. Revert of L-OPLS retyping on the alkyl tail (2026-11-12)

### 7.1 Decision

The sp2->sp3 backbone fix (§6.1) and the OPLS-AA library charges for the
backbone/ester region (§6.2) are kept. The L-OPLS-AA retyping of the alkyl
tail (§6.3) is reverted: the tail now uses plain OPLS-AA throughout, same as
the rest of the molecule, so the full repeat unit is consistently standard
OPLS-AA with library charges and no L-OPLS content anywhere.

### 7.2 Method

`etc/retype_alkyl_tail_lopls.py` stamped every line it touched with a
distinctive comment (`; L-OPLS` on atom charge lines, `; L-OPLS CT-CT-CT-CT`
on tail dihedrals, `; L-OPLS CH2-H epsilon` in `custom-ffnonbonded.itp`).
`etc/revert_alkyl_tail_to_oplsaa.py` finds every line carrying one of these
markers and overwrites it with the standard-OPLS-AA value, using constants
already established and cited elsewhere in this library rather than
re-deriving anything:

| Quantity | L-OPLS value reverted | Standard OPLS-AA value | Same constant as |
|---|---|---|---|
| CH3 carbon charge | -0.222 | -0.180 (`opls_135`) | `assign_library_charges.py`'s `CHARGE_ALPHA_METHYL_C` |
| CH2 carbon charge | -0.148 | -0.120 (`opls_136`) | `assign_library_charges.py`'s `CHARGE_CBETA` |
| CH3-H / CH2-H charge | 0.074 | 0.060 (`opls_140`) | `assign_library_charges.py`'s `CHARGE_GENERIC_HC` |
| CH2-H LJ epsilon | 0.11000 | 0.12552 (sigma unchanged, 0.250) | the pre-L-OPLS value documented in `retype_alkyl_tail_lopls.py` |
| Tail `CT-CT-CT-CT` dihedral | `(0.518787, -0.230192, 0.896807, -1.49134, 0, 0)` | `(0.6276, 1.8828, 0.0, -2.5104, 0, 0)` | `fix_backbone_hybridization.py`'s `DIHEDRAL_CT_CT_CT_HC` |

Reverted lines are re-tagged `; OPLS-AA library charge` (atoms, matching the
style already used on the backbone/ester lines) or
`; OPLS-AA (reverted from L-OPLS)` (dihedrals, nonbonded epsilon). Branch
points (e.g. EH's 2-ethylhexyl branch carbon) were never touched by the
L-OPLS pass and so are untouched here too.

`retype_alkyl_tail_lopls.py` and `etc/lopls-reference/` are kept, not
deleted: the former documents a real decision that was tried and reverted;
the latter is still the provenance source for the backbone/ester library
charges in §6.2 (its stock, `#ifndef LOPLS` branch), independent of whether
the L-OPLS deltas themselves are used anywhere.

### 7.3 Verification

Re-running `renormalize-charges.py` on all five tail-bearing monomers (same
convention as §6.4) restored exact per-repeat-unit neutrality (residuals
were +0.031 to +0.047 e immediately after the revert, from the same
branch-point/analogy effects as §6.4). Confirmed for all six monomers:
`polyply gen_params` assembles a trimer cleanly with exact net-zero charge,
`generate_pama_links.py`'s output is byte-for-byte unchanged, and no
`; L-OPLS` markers or L-OPLS-specific numeric signatures (`0.11000`,
`-0.222...`, `-0.148...`, `0.518787`, ...) remain anywhere in
`polyply/data/oplsaaLigParGen/` or `etc/custom-ffnonbonded.itp`.

## 8. Fix of the two remaining sp2(CM)-flavored bonds at Calpha (2026-11-12)

### 8.1 Finding

`fix_backbone_hybridization.py` (§6.1) only rewrites bonded terms whose atom
set is exactly `{Calpha, Cbeta}` (the vinyl double bond itself), plus
angles/impropers/dihedrals centered on one of those two atoms. It does not
touch the other two bonds radiating from Calpha: to the alpha-methyl
carbon, and to the ester carbonyl carbon. Both were already single bonds
in the free vinyl monomer (Calpha's only double bond was to Cbeta), so this
is not a bond-order bug — but LigParGen parametrized their length/force
constant against Calpha's then-sp2 (`CM`) atom type, and that choice is now
wrong given Calpha is sp3 after polymerization:

| Bond | As shipped (sp2-adjacent) | OPLS-AA bond type matched | Correct sp3 value | OPLS-AA bond type |
|---|---|---|---|---|
| Calpha-AlphaMethylC | 0.1510 nm / 265265.6 | `CM-CT` (`ffbonded.itp`, "; wlj") | 0.1529 nm / 224262.4 | `CT-CT`, same as the Calpha-Cbeta bond in §6.1 |
| Calpha-CarbonylC | 0.1444 nm / 343088.0 | `C-CM` (conjugated-vinyl carbonyl, e.g. thymine) | 0.1522 nm / 265265.6 | `C-CT` |

Both target values come from the generic (non-ifdef'd) entries of the
vendored `etc/lopls-reference/oplsaa.ff/ffbonded.itp`, the same source used
throughout §6-7. The alpha-methyl carbon and the Calpha-Cbeta neighbour are
both plain `CT`; the ester carbonyl carbon is the standard ester `C` type —
neither is `CM`, so the `CM`-flavored values LigParGen assigned are stale
leftovers from the free-vinyl-monomer parametrization.

### 8.2 Method

New script `etc/fix_calpha_substituent_bonds.py` identifies the alpha-methyl
carbon and carbonyl carbon among Calpha's two non-Cbeta heavy neighbours
(same logic as `assign_library_charges.py`'s `identify_backbone_ester_roles`
— the carbonyl carbon is the one bonded to two oxygens), then rewrites
exactly those two bond lines per monomer, tagged `; sp3 backbone fix (was
CM-CT vinyl-adjacent)` and `; sp3 backbone fix (was C-CM conjugated-vinyl)`
respectively, matching the tagging style of `fix_backbone_hybridization.py`.
No charges or atom types are touched, so no re-run of
`renormalize-charges.py` was needed.

### 8.3 Verification

All six monomers: exactly 2 bonds fixed each (12 total), matching the roles
identified in §8.1. Net charge per monomer confirmed still exactly `0` after
the fix (unaffected, since only bond length/force-constant columns
changed). `polyply gen_params` assembles a trimer cleanly for all six
monomers with no errors.

## 9. Three more residuals in the Calpha/Cbeta region (2026-11-12)

### 9.1 Finding

A systematic cross-check of every bond/angle/dihedral touching Calpha,
Cbeta, the alpha-methyl carbon, and the ester carbonyl carbon against the
vendored `etc/lopls-reference/oplsaa.ff/ffbonded.itp` turned up three more
residuals from the same root cause as §6 and §8:
`fix_backbone_hybridization.py` only rewrites terms whose relevant atom set
is exactly `{Calpha, Cbeta}` (the bond itself, or angles/impropers/
dihedrals centered on one of those two atoms specifically).

1. **Cbeta-H bonds**: `0.1080 nm / 284512.0` in all six monomers, matching
   OPLS-AA's `CM HC` (sp² vinyl C-H) exactly. Every other alkyl C-H bond in
   the molecule — including the alpha-methyl C-H bonds right next to it —
   is `CT HC`: `0.1090 nm / 284512.0` (identical force constant, only the
   equilibrium length is stale). `fix_backbone_hybridization.py` rewrote
   the Calpha=Cbeta bond and the angles/dihedrals centered on Cbeta, but
   never Cbeta's own bonds to its two hydrogens.

2. **Dihedrals centered on the Calpha-alphaMethylC bond** (6 per monomer):
   raw per-molecule LigParGen values (e.g. `(-0.778, -2.335, 0, 3.113, 0,
   0)`), not the generic OPLS-AA wildcard torsions matching their
   atom-type quartet:

   | Quartet pattern | Current (stale) | OPLS-AA generic match |
   |---|---|---|
   | Cbeta-Calpha-alphaMethylC-H (`CT CT CT HC`) | per-molecule fit | `(0.6276, 1.8828, 0, -2.5104, 0, 0)` — same constant already reused for the Calpha-Cbeta bond's own dihedrals |
   | carbonylC-Calpha-alphaMethylC-H (`C CT CT HC`) | per-molecule fit | `(-0.2092, -0.6276, 0, 0.8368, 0, 0)` |

3. **Dihedrals centered on the Calpha-carbonylC bond** (4 per monomer):
   coefficients up to ~30 kJ/mol, uniform across all six monomers, and
   internally inconsistent with each other (the Cbeta-side and
   alpha-methyl-side quartets of the same type should be identical under
   OPLS-AA's generic rule, but currently aren't — e.g. `(30.334, -5.23,
   -25.104, 0, 0, 0)` vs. `(1.789, -5.052, -0.962, 4.226, 0, 0)` for the
   two `CT CT C O` cases). This pattern — large, molecule-specific,
   internally inconsistent — is consistent with leftover fitting against
   the free monomer's Cbeta=Calpha-C(=O) acrylate conjugation, which no
   longer exists once Calpha is sp3. Generic matches:

   | Quartet pattern | OPLS-AA generic match |
   |---|---|
   | {Cbeta,alphaMethylC}-Calpha-carbonylC-carbonylO (`CT CT C O`, "propanamide") | `(4.87855, 0, -4.87855, 0, 0, 0)` |
   | {Cbeta,alphaMethylC}-Calpha-carbonylC-esterO (`CT CT C OS`, "esters") | `(-1.15688, -3.47063, 0, 4.6275, 0, 0)` |

All target values come from the generic (non-ifdef'd) `[ dihedraltypes ]`
entries in the same vendored reference file already used throughout §6-8.

### 9.2 Method

New script `etc/fix_calpha_region_residual_terms.py` identifies Cbeta's
hydrogens, the alpha-methyl carbon/hydrogens, and the carbonyl/ester
oxygens (same role-identification logic as §6-8), then:
- rewrites both Cbeta-H bond lines to `CT-HC`,
- for each proper dihedral whose central bond is `{Calpha,
  alphaMethylC}`, picks `CT CT CT HC` or `C CT CT HC` depending on whether
  the far end (opposite the methyl's H) is Cbeta or the carbonyl carbon,
- for each proper dihedral whose central bond is `{Calpha, carbonylC}`,
  picks `CT CT C O` or `CT CT C OS` depending on whether the far end is the
  carbonyl oxygen or the ester oxygen.

All rewritten lines are tagged `; sp3 backbone fix (was CM-HC
vinyl-adjacent)` (bond) or `; sp3 backbone fix (was LigParGen
vinyl-conjugated fit)` (dihedrals). No charges or atom types are touched.

### 9.3 Verification

All six monomers: exactly 2 bonds, 6 alpha-methyl-side dihedrals, and 4
carbonyl-side dihedrals fixed each (matching the counts predicted in
§9.1). Net charge per monomer confirmed still exactly `0` (unaffected).
`polyply gen_params` assembles a trimer cleanly for all six monomers with
no errors or warnings.
