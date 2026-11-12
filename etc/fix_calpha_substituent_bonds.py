#!/usr/bin/env python3
"""Fix the two remaining sp2(CM)-flavored bonds radiating from Calpha.

fix_backbone_hybridization.py (2026-10-28) rewrote the Calpha=Cbeta vinyl
double bond itself, plus the angles/impropers/dihedrals centered on
Calpha/Cbeta, to sp3 values -- but it only touched the bond whose atom pair
was exactly {Calpha, Cbeta}. Two other bonds radiating from Calpha were left
untouched. They were already single bonds in the free vinyl monomer too, so
this isn't a bond-order fix, but their equilibrium length/force constant
were parametrized by LigParGen against Calpha's then-sp2 (CM) identity:

  - Calpha-AlphaMethylC: 0.1510 nm / 265265.6 kJ/mol/nm^2, matching OPLS-AA's
    "CM CT" bond type (vendored etc/lopls-reference/oplsaa.ff/ffbonded.itp,
    "; wlj"). Now that Calpha is sp3 (plain CT, same as Cbeta and the
    alpha-methyl carbon), this should be the generic "CT CT" value already
    used for the Calpha-Cbeta bond: 0.1529 nm / 224262.4.
  - Calpha-CarbonylC: 0.1444 nm / 343088.0, matching "C CM" (a
    conjugated-vinyl carbonyl bond, e.g. thymine's C5=C6-C4=O system, in the
    same vendored file). The correct value for an sp3 carbon bonded to an
    ester carbonyl carbon is "C CT": 0.1522 nm / 265265.6.

Role identification (alpha-methyl C vs. carbonyl C among Calpha's two
non-Cbeta heavy neighbours) is the same logic as
etc/assign_library_charges.py's identify_backbone_ester_roles.

Usage:
    python fix_calpha_substituent_bonds.py
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

BOND_CALPHA_ALPHAMETHYL = (0.15290, 224262.400)  # CT-CT
BOND_CALPHA_CARBONYL = (0.15220, 265265.600)  # C-CT


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


def identify_alpha_methyl_and_carbonyl(idx_to_name, adjacency, calpha, cbeta):
    heavy = [n for n in adjacency[calpha] if not is_h(idx_to_name, n)]
    others = [n for n in heavy if n != cbeta]
    assert len(others) == 2

    def oxygen_count(idx):
        return sum(1 for m in adjacency[idx] if idx_to_name[m].startswith("O"))

    carbonyl_c, alpha_methyl_c = (
        others if oxygen_count(others[0]) > oxygen_count(others[1]) else others[::-1]
    )
    return alpha_methyl_c, carbonyl_c


def fix_monomer(path, calpha_name, cbeta_name):
    lines = path.read_text().splitlines(keepends=True)
    idx_to_name, name_to_idx, adjacency = parse_structure([l.rstrip("\n") for l in lines])
    calpha = name_to_idx[calpha_name]
    cbeta = name_to_idx[cbeta_name]
    alpha_methyl_c, carbonyl_c = identify_alpha_methyl_and_carbonyl(
        idx_to_name, adjacency, calpha, cbeta
    )

    out = []
    n_fixed = 0
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
            if pair == {calpha, alpha_methyl_c}:
                out.append(
                    "{:6d}{:6d}{:6d}{:12.4f}{:12.3f}  ; sp3 backbone fix (was CM-CT vinyl-adjacent)\n".format(
                        ai, aj, 1, *BOND_CALPHA_ALPHAMETHYL
                    )
                )
                n_fixed += 1
                continue
            if pair == {calpha, carbonyl_c}:
                out.append(
                    "{:6d}{:6d}{:6d}{:12.4f}{:12.3f}  ; sp3 backbone fix (was C-CM conjugated-vinyl)\n".format(
                        ai, aj, 1, *BOND_CALPHA_CARBONYL
                    )
                )
                n_fixed += 1
                continue

        out.append(line)

    path.write_text("".join(out))
    print(
        f"{path.name}: bonds_fixed={n_fixed} "
        f"(Calpha={calpha_name}, AlphaMethylC={idx_to_name[alpha_methyl_c]}, "
        f"CarbonylC={idx_to_name[carbonyl_c]})"
    )


def main():
    for mono, (calpha_name, cbeta_name) in MONOMERS.items():
        path = LIB_DIR / f"PAMA_{mono}.oplsaa.LigParGen.ff"
        fix_monomer(path, calpha_name, cbeta_name)


if __name__ == "__main__":
    main()
