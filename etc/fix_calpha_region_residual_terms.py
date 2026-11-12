#!/usr/bin/env python3
"""Fix the remaining sp2(vinyl)-flavored terms in the Calpha/Cbeta region:
the Cbeta-H bond length, and the proper dihedrals centered on the two
Calpha bonds that fix_backbone_hybridization.py / fix_calpha_substituent_
bonds.py already fixed at the *bond* level but not at the *dihedral* level.

fix_backbone_hybridization.py only rewrites dihedrals whose central bond
is exactly {Calpha, Cbeta}. It also never touched Cbeta's own bonds to its
two hydrogens. Three residuals are left over from the free vinyl monomer's
sp2/conjugated topology, found by cross-matching current values against
the vendored etc/lopls-reference/oplsaa.ff/ffbonded.itp:

1. Cbeta-H bonds: 0.1080 nm / 284512.0, matching OPLS-AA's `CM HC` (vinyl
   =CH2 C-H). Every other alkyl C-H bond in the molecule, including the
   alpha-methyl C-H bonds right next to it, is `CT HC`: 0.1090 nm /
   284512.0 (same force constant, only the equilibrium length is stale).

2. Dihedrals with central bond Calpha-AlphaMethylC: raw per-molecule
   LigParGen fits instead of the generic wildcard torsions that match
   their atom-type quartet:
     - {Cbeta,AlphaMethylC-H}  -> `CT CT CT HC` = (0.6276, 1.8828, 0, -2.5104, 0, 0)
       (the same constant already reused for the Calpha-Cbeta bond itself)
     - {CarbonylC,AlphaMethylC-H} -> `C CT CT HC` = (-0.2092, -0.6276, 0, 0.8368, 0, 0)

3. Dihedrals with central bond Calpha-CarbonylC: large (~30 kJ/mol)
   per-molecule coefficients consistent with the now-obsolete
   Cbeta=Calpha-C(=O) acrylate conjugation (and inconsistent with each
   other across the Cbeta-side vs. alpha-methyl-side quartets, which OPLS-
   AA's generic rule says should be identical). Replaced with:
     - {Cbeta or AlphaMethylC, CarbonylO} -> `CT CT C O` ("propanamide")
       = (4.87855, 0, -4.87855, 0, 0, 0)
     - {Cbeta or AlphaMethylC, EsterO}    -> `CT CT C OS` ("esters")
       = (-1.15688, -3.47063, 0, 4.6275, 0, 0)

See etc/report.md section 9 for the full rationale and verification.

Usage:
    python fix_calpha_region_residual_terms.py
"""
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
LIB_DIR = REPO_ROOT / "polyply" / "data" / "oplsaaLigParGen"

MONOMERS = {
    "EH": ("C04", "C05"),
    "LA": ("C04", "C05"),
    "MMA": ("C07", "C08"),
    "OC": ("C03", "C01"),
    "ST": ("C04", "C05"),
    "TD": ("C18", "C19"),
}

BOND_CBETA_H = (0.10900, 284512.000)  # CT-HC

DIH_CT_CT_CT_HC = (0.627600, 1.882800, 0.000000, -2.510400, 0.000000, 0.000000)
DIH_C_CT_CT_HC = (-0.209200, -0.627600, 0.000000, 0.836800, 0.000000, 0.000000)
DIH_CT_CT_C_O = (4.878550, 0.000000, -4.878550, 0.000000, 0.000000, 0.000000)
DIH_CT_CT_C_OS = (-1.156880, -3.470630, 0.000000, 4.627500, 0.000000, 0.000000)


def parse_structure(lines):
    idx_to_name, name_to_idx, adjacency = {}, {}, {}
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
            nr, name = int(t[0]), t[4]
            idx_to_name[nr] = name
            name_to_idx[name] = nr
            adjacency[nr] = set()
        elif in_bonds:
            t = s.split()
            a, b = int(t[0]), int(t[1])
            adjacency.setdefault(a, set()).add(b)
            adjacency.setdefault(b, set()).add(a)
    return idx_to_name, name_to_idx, adjacency


def is_h(idx_to_name, idx):
    return idx_to_name[idx].startswith("H")


def identify_roles(idx_to_name, adjacency, calpha, cbeta):
    cbeta_h = [n for n in adjacency[cbeta] if is_h(idx_to_name, n)]
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
    return {
        "calpha": calpha, "cbeta": cbeta, "cbeta_h": cbeta_h,
        "carbonyl_c": carbonyl_c, "alpha_methyl_c": alpha_methyl_c,
        "alpha_methyl_h": alpha_methyl_h,
        "carbonyl_o": carbonyl_o, "ester_o": ester_o,
    }


def fix_monomer(path, calpha_name, cbeta_name):
    lines = path.read_text().splitlines(keepends=True)
    idx_to_name, name_to_idx, adjacency = parse_structure([l.rstrip("\n") for l in lines])
    calpha = name_to_idx[calpha_name]
    cbeta = name_to_idx[cbeta_name]
    roles = identify_roles(idx_to_name, adjacency, calpha, cbeta)

    out = []
    n_bond = n_dih_methyl = n_dih_carbonyl = 0
    atoms_section_done = False
    for line in lines:
        stripped = line.strip()
        if not atoms_section_done:
            out.append(line)
            if stripped.startswith("[ bonds ]"):
                atoms_section_done = True
            continue

        tokens = stripped.split()
        is_data_line = (
            bool(tokens)
            and not stripped.startswith(";")
            and not stripped.startswith("[")
            and all(re.match(r"^-?[0-9.eE+-]+$", t) for t in tokens[: min(4, len(tokens))])
        )

        if is_data_line and len(tokens) == 5 and int(tokens[2]) == 1:
            ai, aj = int(tokens[0]), int(tokens[1])
            pair = {ai, aj}
            if roles["cbeta"] in pair:
                other = aj if ai == roles["cbeta"] else ai
                if other in roles["cbeta_h"]:
                    out.append(
                        "{:6d}{:6d}{:6d}{:12.4f}{:12.3f}  ; sp3 backbone fix (was CM-HC vinyl-adjacent)\n".format(
                            ai, aj, 1, *BOND_CBETA_H
                        )
                    )
                    n_bond += 1
                    continue

        if is_data_line and len(tokens) in (10, 11) and int(tokens[4]) == 3:
            ai, aj, ak, al = int(tokens[0]), int(tokens[1]), int(tokens[2]), int(tokens[3])
            central = {aj, ak}

            if central == {roles["calpha"], roles["alpha_methyl_c"]}:
                ends = [ai, al]
                h_end = next((e for e in ends if is_h(idx_to_name, e)), None)
                heavy_end = next((e for e in ends if e != h_end), None)
                if h_end is not None and heavy_end is not None:
                    target = DIH_CT_CT_CT_HC if heavy_end == roles["cbeta"] else DIH_C_CT_CT_HC
                    out.append(
                        "{:6d}{:6d}{:6d}{:6d}{:6d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}"
                        "  ; sp3 backbone fix (was LigParGen vinyl-conjugated fit)\n".format(
                            ai, aj, ak, al, 3, *target
                        )
                    )
                    n_dih_methyl += 1
                    continue

            if central == {roles["calpha"], roles["carbonyl_c"]}:
                ends = [ai, al]
                o_end = next((e for e in ends if e in (roles["carbonyl_o"], roles["ester_o"])), None)
                if o_end is not None:
                    target = DIH_CT_CT_C_O if o_end == roles["carbonyl_o"] else DIH_CT_CT_C_OS
                    out.append(
                        "{:6d}{:6d}{:6d}{:6d}{:6d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}"
                        "  ; sp3 backbone fix (was LigParGen vinyl-conjugated fit)\n".format(
                            ai, aj, ak, al, 3, *target
                        )
                    )
                    n_dih_carbonyl += 1
                    continue

        out.append(line)

    path.write_text("".join(out))
    print(
        f"{path.name}: cbeta_h_bonds_fixed={n_bond} "
        f"alphaMethylC_dihedrals_fixed={n_dih_methyl} carbonylC_dihedrals_fixed={n_dih_carbonyl}"
    )


def main():
    for mono, (calpha_name, cbeta_name) in MONOMERS.items():
        path = LIB_DIR / f"PAMA_{mono}.oplsaa.LigParGen.ff"
        fix_monomer(path, calpha_name, cbeta_name)


if __name__ == "__main__":
    main()
