#!/usr/bin/env python3
"""Fix the remaining OPLS-AA residuals found by a full bond/angle/dihedral/
LJ audit of every PAMA monomer against the vendored stock OPLS-AA tables
(etc/lopls-reference/oplsaa.ff/{ffbonded,ffnonbonded}.itp).

Six issues, verified by direct grep of the reference tables:

1. LJ: the ester carbonyl carbon (role "C") keeps LigParGen's raw CM1A fit
   (sigma=0.355, epsilon=0.29288) in etc/custom-ffnonbonded.itp, even though
   assign_library_charges.py already switched its *charge* to the opls_465
   value. Correct LJ (opls_465): sigma=0.375, epsilon=0.43932.

2. LJ: MMA's 3 methoxy hydrogens keep generic HC LJ (sigma=0.25,
   epsilon=0.12552) even though their charge was already switched to the
   opls_469 value. Correct LJ (opls_469): sigma=0.242, epsilon=0.06276.
   (The other monomers' first-alkoxy-C is a CH2, not a methoxy CH3, and has
   no dedicated H LJ entry in stock OPLS-AA -- generic HC there is already
   the documented best-available analogy, not a bug; see report.md 6.2.)

3. Angle: alpha-methyl H - alpha-methyl C - Calpha (HC-CT-CT) is still a raw
   LigParGen fit (109.500/292.880). Correct generic value: 110.700/313.800.

4. Angle: fix_backbone_hybridization.py's original blanket angle fix set
   *every* angle centered on Calpha/Cbeta to the symmetric CT-CT-CT value
   (112.700/488.273), including the two where the outer atom is really the
   ester carbonyl carbon (role "C", not "CT"): {Cbeta,AlphaMethylC}-Calpha-
   CarbonylC. Those need the distinct CT-CT-C value: 111.100/527.184.

5. Angle: Calpha-CarbonylC-CarbonylO (CT-C-O) is still a raw LigParGen fit
   at 125.300 (same k=669.440 as the correct generic entry, 120.400/
   669.440) -- a residual trace of the free monomer's acrylate conjugation
   on the carbonyl's own geometry.

6. Dihedral: the same blanket-fix issue as (4), one bond over. The original
   fix_backbone_hybridization.py dihedral fix applied the CT-CT-CT-HC
   wildcard to all four dihedrals around the Calpha-Cbeta bond, including
   the two where the far atom is CarbonylC (role "C"). Those need
   C-CT-CT-HC = (-0.2092,-0.6276,0,0.8368,0,0) -- the same correction
   already applied to the sibling Calpha-AlphaMethylC bond by
   fix_calpha_region_residual_terms.py.

7. Dihedral: revert_alkyl_tail_to_oplsaa.py's tail CT-CT-CT-CT dihedral
   (four consecutive backbone carbons, no H at either end) was given the
   wrong target: it reused (0.6276,1.8828,0,-2.5104,0,0), which is actually
   OPLS-AA's CT-CT-CT-HC torsion (H at one end), not CT-CT-CT-CT. The
   correct all-carbon value is (2.9288,-1.4644,0.2092,-1.6736,0,0) -- a
   different curve shape, not just different coefficients. Identified by
   its existing "; OPLS-AA (reverted from L-OPLS)" tag.

8. Dihedral: the ester's own CT-OS-C-O torsion (first-alkoxy-C - esterO -
   carbonylC - carbonylO) is a raw two-term LigParGen fit, (21.439,0,
   -21.439,0,0,0), close to but not equal to the one stock entry for this
   exact role-tuple, (20.92,0,-20.92,0,0,0; "benzoic acids & esters").

See etc/report.md section 10 for the full rationale and verification.

Usage:
    python fix_remaining_opls_residuals.py
"""
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
LIB_DIR = REPO_ROOT / "polyply" / "data" / "oplsaaLigParGen"
NB_FILE = REPO_ROOT / "etc" / "custom-ffnonbonded.itp"

MONOMERS = {
    "EH": ("C04", "C05"),
    "LA": ("C04", "C05"),
    "MMA": ("C07", "C08"),
    "OC": ("C03", "C01"),
    "ST": ("C04", "C05"),
    "TD": ("C18", "C19"),
}

LJ_CARBONYL_C = (0.37500, 0.43932)      # opls_465
LJ_METHOXY_H = (0.24200, 0.06276)       # opls_469

ANGLE_HC_CT_CT = (110.700, 313.800)
ANGLE_CT_CT_C = (111.100, 527.184)
ANGLE_CT_C_O = (120.400, 669.440)

DIH_C_CT_CT_HC = (-0.209200, -0.627600, 0.000000, 0.836800, 0.000000, 0.000000)
DIH_CT_CT_CT_CT = (2.928800, -1.464400, 0.209200, -1.673600, 0.000000, 0.000000)
DIH_CT_OS_C_O = (20.920000, 0.000000, -20.920000, 0.000000, 0.000000, 0.000000)


def parse_structure(lines):
    idx_to_name, name_to_idx, idx_to_type, adjacency = {}, {}, {}, {}
    in_atoms, in_bonds = False, False
    for line in lines:
        s = line.strip()
        if s.startswith("[ atoms ]"):
            in_atoms, in_bonds = True, False
            continue
        if s.startswith("[ bonds ]"):
            in_atoms, in_bonds = False, True
            continue
        if s.startswith("["):
            if in_bonds:
                break
            in_atoms = False
            continue
        if not s or s.startswith(";"):
            continue
        if in_atoms:
            t = s.split()
            nr, name, atype = int(t[0]), t[4], t[1]
            idx_to_name[nr] = name
            name_to_idx[name] = nr
            idx_to_type[nr] = atype
            adjacency[nr] = set()
        elif in_bonds:
            t = s.split()
            a, b = int(t[0]), int(t[1])
            adjacency.setdefault(a, set()).add(b)
            adjacency.setdefault(b, set()).add(a)
    return idx_to_name, name_to_idx, idx_to_type, adjacency


def is_h(idx_to_name, idx):
    return idx_to_name[idx].startswith("H")


def identify_roles(idx_to_name, adjacency, calpha, cbeta):
    heavy = [n for n in adjacency[calpha] if not is_h(idx_to_name, n)]
    others = [n for n in heavy if n != cbeta]

    def oxygen_count(idx):
        return sum(1 for m in adjacency[idx] if idx_to_name[m].startswith("O"))

    carbonyl_c, alpha_methyl_c = (
        others if oxygen_count(others[0]) > oxygen_count(others[1]) else others[::-1]
    )
    alpha_methyl_h = [n for n in adjacency[alpha_methyl_c] if is_h(idx_to_name, n)]
    carbonyl_oxygens = [n for n in adjacency[carbonyl_c] if idx_to_name[n].startswith("O")]
    ester_o, carbonyl_o = (
        carbonyl_oxygens
        if len([m for m in adjacency[carbonyl_oxygens[0]] if not is_h(idx_to_name, m)]) == 2
        else carbonyl_oxygens[::-1]
    )
    first_alkoxy_c = next(n for n in adjacency[ester_o] if n != carbonyl_c)
    first_alkoxy_h = [n for n in adjacency[first_alkoxy_c] if is_h(idx_to_name, n)]
    return {
        "calpha": calpha, "cbeta": cbeta,
        "carbonyl_c": carbonyl_c, "alpha_methyl_c": alpha_methyl_c,
        "alpha_methyl_h": alpha_methyl_h,
        "carbonyl_o": carbonyl_o, "ester_o": ester_o,
        "first_alkoxy_c": first_alkoxy_c, "first_alkoxy_h": first_alkoxy_h,
    }


def is_data_line(tokens, stripped, min_numeric_cols):
    return (
        bool(tokens)
        and not stripped.startswith(";")
        and not stripped.startswith("[")
        and all(re.match(r"^-?[0-9.eE+-]+$", t) for t in tokens[:min_numeric_cols])
    )


def fix_monomer(path, calpha_name, cbeta_name):
    lines = path.read_text().splitlines(keepends=True)
    idx_to_name, name_to_idx, idx_to_type, adjacency = parse_structure(
        [l.rstrip("\n") for l in lines]
    )
    calpha = name_to_idx[calpha_name]
    cbeta = name_to_idx[cbeta_name]
    roles = identify_roles(idx_to_name, adjacency, calpha, cbeta)

    carbonyl_c_type = idx_to_type[roles["carbonyl_c"]]
    methoxy_h_types = set()
    if len(roles["first_alkoxy_h"]) == 3:  # terminal CH3 methoxy (MMA only)
        methoxy_h_types = {idx_to_type[h] for h in roles["first_alkoxy_h"]}

    out = []
    n_angle_hc_ct_ct = n_angle_ct_ct_c = n_angle_ct_c_o = 0
    n_dih_carbonyl_cbeta = n_dih_tail = n_dih_ester = 0
    atoms_section_done = False
    for line in lines:
        stripped = line.strip()
        if not atoms_section_done:
            out.append(line)
            if stripped.startswith("[ bonds ]"):
                atoms_section_done = True
            continue

        tokens = stripped.split(";", 1)[0].split()

        # Angles: 6 cols (ai aj ak funct c0 c1), funct==1
        if is_data_line(tokens, stripped, 4) and len(tokens) == 6 and tokens[3] == "1":
            ai, aj, ak = int(tokens[0]), int(tokens[1]), int(tokens[2])
            triple = {ai, ak}

            # (3) alpha-methyl H - alpha-methyl C - Calpha, raw LigParGen value
            if aj == roles["alpha_methyl_c"] and roles["calpha"] in triple:
                h_end = ai if ai in roles["alpha_methyl_h"] else (ak if ak in roles["alpha_methyl_h"] else None)
                if h_end is not None:
                    out.append(
                        "{:6d}{:6d}{:6d}{:6d}{:12.3f}{:12.3f}  ; corrected HC-CT-CT (was raw LigParGen fit)\n".format(
                            ai, aj, ak, 1, *ANGLE_HC_CT_CT
                        )
                    )
                    n_angle_hc_ct_ct += 1
                    continue

            # (4) {Cbeta,AlphaMethylC}-Calpha-CarbonylC, blanket-fixed to CT-CT-CT
            if aj == roles["calpha"] and roles["carbonyl_c"] in triple:
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:12.3f}{:12.3f}  ; corrected CT-CT-C (was blanket CT-CT-CT backbone fix)\n".format(
                        ai, aj, ak, 1, *ANGLE_CT_CT_C
                    )
                )
                n_angle_ct_ct_c += 1
                continue

            # (5) Calpha-CarbonylC-CarbonylO, raw LigParGen value
            if aj == roles["carbonyl_c"] and roles["calpha"] in triple and roles["carbonyl_o"] in triple:
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:12.3f}{:12.3f}  ; corrected CT-C-O (was raw LigParGen fit)\n".format(
                        ai, aj, ak, 1, *ANGLE_CT_C_O
                    )
                )
                n_angle_ct_c_o += 1
                continue

        # Dihedrals: 11 cols (ai aj ak al funct c0..c5), funct==3
        if is_data_line(tokens, stripped, 4) and len(tokens) == 11 and tokens[4] == "3":
            ai, aj, ak, al = int(tokens[0]), int(tokens[1]), int(tokens[2]), int(tokens[3])
            central = {aj, ak}
            ends = {ai, al}

            # (6) CarbonylC-Calpha-Cbeta-H, blanket-fixed to CT-CT-CT-HC
            if central == {roles["calpha"], roles["cbeta"]} and roles["carbonyl_c"] in ends:
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:6d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}"
                    "  ; corrected C-CT-CT-HC (was blanket CT-CT-CT-HC backbone fix)\n".format(
                        ai, aj, ak, al, 3, *DIH_C_CT_CT_HC
                    )
                )
                n_dih_carbonyl_cbeta += 1
                continue

            # (7) tail CT-CT-CT-CT, mis-set to the CT-CT-CT-HC value
            if stripped.rstrip().endswith("; OPLS-AA (reverted from L-OPLS)"):
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:6d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}"
                    "  ; OPLS-AA (reverted from L-OPLS, corrected CT-CT-CT-CT value)\n".format(
                        ai, aj, ak, al, 3, *DIH_CT_CT_CT_CT
                    )
                )
                n_dih_tail += 1
                continue

            # (8) ester CT-OS-C-O torsion, raw LigParGen value
            if (
                {roles["first_alkoxy_c"], roles["carbonyl_o"]} == ends
                and roles["ester_o"] in central
                and roles["carbonyl_c"] in central
            ):
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:6d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}"
                    "  ; corrected CT-OS-C-O (was raw LigParGen fit)\n".format(
                        ai, aj, ak, al, 3, *DIH_CT_OS_C_O
                    )
                )
                n_dih_ester += 1
                continue

        out.append(line)

    path.write_text("".join(out))
    print(
        f"{path.name}: HC-CT-CT_angles={n_angle_hc_ct_ct} CT-CT-C_angles={n_angle_ct_ct_c} "
        f"CT-C-O_angles={n_angle_ct_c_o} C-CT-CT-HC_dihedrals={n_dih_carbonyl_cbeta} "
        f"tail_CT-CT-CT-CT_dihedrals={n_dih_tail} ester_dihedrals={n_dih_ester}"
    )
    return carbonyl_c_type, methoxy_h_types


def fix_nonbonded(carbonyl_c_types, methoxy_h_types):
    text = NB_FILE.read_text()
    out, n_carbonyl, n_methoxy = [], 0, 0
    for line in text.splitlines(keepends=True):
        t = line.split()
        if t and t[0] in carbonyl_c_types:
            name, mass, charge, ptype = t[0], float(t[1]), t[2], t[3]
            out.append(
                "{:<12s}{:>10.4f}    {:>5s}    {} {:.5E} {:.5E}  ; corrected to opls_465 LJ (was raw LigParGen fit)\n".format(
                    name, mass, charge, ptype, *LJ_CARBONYL_C
                )
            )
            n_carbonyl += 1
        elif t and t[0] in methoxy_h_types:
            name, mass, charge, ptype = t[0], float(t[1]), t[2], t[3]
            out.append(
                "{:<12s}{:>10.4f}    {:>5s}    {} {:.5E} {:.5E}  ; corrected to opls_469 LJ (was generic HC)\n".format(
                    name, mass, charge, ptype, *LJ_METHOXY_H
                )
            )
            n_methoxy += 1
        else:
            out.append(line)
    NB_FILE.write_text("".join(out))
    print(f"custom-ffnonbonded.itp: carbonyl_C_LJ_fixed={n_carbonyl} methoxy_H_LJ_fixed={n_methoxy}")


def main():
    carbonyl_c_types = set()
    methoxy_h_types = set()
    for mono, (calpha_name, cbeta_name) in MONOMERS.items():
        path = LIB_DIR / f"PAMA_{mono}.oplsaa.LigParGen.ff"
        c_type, h_types = fix_monomer(path, calpha_name, cbeta_name)
        carbonyl_c_types.add(c_type)
        methoxy_h_types |= h_types
    fix_nonbonded(carbonyl_c_types, methoxy_h_types)


if __name__ == "__main__":
    main()
