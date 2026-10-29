#!/usr/bin/env python3
"""Fix the sp2 (vinyl C=C) backbone bug in the PAMA monomer library.

LigParGen parametrized each PAMA monomer as the free vinyl monomer, so the
intra-residue Calpha=Cbeta bond is still a double bond (OPLS CM=CM values),
with sp2-planarity improper dihedrals at both centres and C=C no-rotation
torsions around that bond. After polymerization this bond is single and both
carbons are sp3 (Calpha is an all-carbon quaternary centre; Cbeta is a plain
backbone CH2). This script rewrites the bonded terms and nonbonded atom types
for Calpha/Cbeta to the generic sp3 OPLS-AA values already used by
generate_pama_links.py for the *inter*-residue link, so intra- and
inter-residue backbone parameters stay consistent.

Usage:
    python fix_backbone_hybridization.py
(no arguments; operates on the fixed monomer table below, run from anywhere)
"""
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
LIB_DIR = REPO_ROOT / "polyply" / "data" / "oplsaaLigParGen"
NONBONDED_FILE = REPO_ROOT / "etc" / "custom-ffnonbonded.itp"

# monomer -> (Calpha atom name, Cbeta atom name), per etc/README.md
MONOMERS = {
    "EH": ("C04", "C05"),
    "LA": ("C04", "C05"),
    "MMA": ("C07", "C08"),
    "OC": ("C03", "C01"),
    "ST": ("C04", "C05"),
    "TD": ("C18", "C19"),
}

# Generic sp3 OPLS-AA values, same source as generate_pama_links.py
BOND_CT_CT = (0.15290, 224262.400)
ANGLE_CT_CT_CT = (112.700, 488.273)
ANGLE_CT_CT_HC = (110.700, 313.800)
ANGLE_HC_CT_HC = (107.800, 276.144)
DIHEDRAL_CT_CT_CT_HC = (0.627600, 1.882800, 0.000000, -2.510400)
NONBONDED_CT = (3.50000e-01, 2.76144e-01)  # sigma, epsilon


def parse_atoms(lines):
    """Return (idx -> name, idx -> type, name -> idx, last_atom_line_idx)."""
    idx_to_name = {}
    idx_to_type = {}
    name_to_idx = {}
    in_atoms = False
    last_idx = 0
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("[ atoms ]"):
            in_atoms = True
            continue
        if in_atoms and stripped.startswith("["):
            break
        if in_atoms:
            if not stripped or stripped.startswith(";"):
                continue
            tokens = stripped.split()
            nr = int(tokens[0])
            atype = tokens[1]
            atom_name = tokens[4]
            idx_to_name[nr] = atom_name
            idx_to_type[nr] = atype
            name_to_idx[atom_name] = nr
            last_idx = i
    return idx_to_name, idx_to_type, name_to_idx, last_idx


def is_h(idx_to_name, idx):
    return idx_to_name[idx].startswith("H")


def fix_monomer(path, calpha_name, cbeta_name):
    lines = path.read_text().splitlines(keepends=True)
    idx_to_name, idx_to_type, name_to_idx, _ = parse_atoms(lines)
    calpha = name_to_idx[calpha_name]
    cbeta = name_to_idx[cbeta_name]
    backbone_pair = {calpha, cbeta}

    out = []
    n_bond = n_angle = n_dih = n_improper_removed = 0
    atoms_section_done = False
    for line in lines:
        stripped = line.strip()
        if not atoms_section_done:
            out.append(line)
            if stripped.startswith("[ bonds ]"):
                atoms_section_done = True
            continue

        tokens = stripped.split()
        is_data_line = bool(tokens) and not stripped.startswith(";") and not stripped.startswith("[") \
            and all(re.match(r"^-?[0-9.eE+-]+$", t) for t in tokens[: min(4, len(tokens))])

        if is_data_line and len(tokens) == 5 and int(tokens[2]) == 1:
            # bonds: ai aj funct b0 kb
            ai, aj = int(tokens[0]), int(tokens[1])
            if {ai, aj} == backbone_pair:
                out.append(
                    "{:6d}{:6d}{:6d}{:12.4f}{:12.3f}  ; sp3 backbone fix (was vinyl C=C)\n".format(
                        ai, aj, 1, *BOND_CT_CT
                    )
                )
                n_bond += 1
                continue

        if is_data_line and len(tokens) == 6 and int(tokens[3]) == 1:
            # angles: ai aj ak funct theta0 k_theta
            ai, aj, ak = int(tokens[0]), int(tokens[1]), int(tokens[2])
            if aj in backbone_pair:
                if aj == calpha:
                    theta0, k = ANGLE_CT_CT_CT
                else:  # aj == cbeta
                    h_count = is_h(idx_to_name, ai) + is_h(idx_to_name, ak)
                    if h_count == 2:
                        theta0, k = ANGLE_HC_CT_HC
                    elif h_count == 1:
                        theta0, k = ANGLE_CT_CT_HC
                    else:
                        theta0, k = ANGLE_CT_CT_CT
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:12.3f}{:12.3f}  ; sp3 backbone fix (was vinyl sp2)\n".format(
                        ai, aj, ak, 1, theta0, k
                    )
                )
                n_angle += 1
                continue

        if is_data_line and len(tokens) == 8 and int(tokens[4]) == 4:
            # improper dihedrals: ai aj ak al funct c0 c1 mult
            aj = int(tokens[1])
            if aj in backbone_pair:
                n_improper_removed += 1
                continue  # drop: neither centre is planar anymore

        if is_data_line and len(tokens) in (10, 11) and int(tokens[4]) == 3:
            # proper RB dihedrals: ai aj ak al funct c0 c1 c2 c3 c4 c5
            aj, ak = int(tokens[1]), int(tokens[2])
            if {aj, ak} == backbone_pair:
                ai, al, funct = int(tokens[0]), int(tokens[3]), int(tokens[4])
                c = list(DIHEDRAL_CT_CT_CT_HC) + [0.000000, 0.000000]
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:6d}{:12.4f}{:12.4f}{:12.4f}{:12.4f}{:12.4f}{:12.4f}"
                    "  ; sp3 backbone fix (was C=C no-rotation)\n".format(
                        ai, aj, ak, al, funct, *c
                    )
                )
                n_dih += 1
                continue

        out.append(line)

    path.write_text("".join(out))
    calpha_type = idx_to_type[calpha]
    cbeta_type = idx_to_type[cbeta]
    print(
        f"{path.name}: bond={n_bond} angles={n_angle} dihedrals={n_dih} "
        f"impropers_removed={n_improper_removed} "
        f"(Calpha={calpha_name}/{calpha_type}, Cbeta={cbeta_name}/{cbeta_type})"
    )
    return calpha_type, cbeta_type


def fix_nonbonded(types_to_fix):
    text = NONBONDED_FILE.read_text()
    lines = text.splitlines(keepends=True)
    out = []
    n_fixed = 0
    for line in lines:
        tokens = line.split()
        if tokens and tokens[0] in types_to_fix:
            name, mass, charge, ptype = tokens[0], float(tokens[1]), tokens[2], tokens[3]
            out.append(
                "{:<12s}{:>10.4f}    {:>5s}    {} {:.5E} {:.5E}  ; sp3 backbone fix (was CM/sp2)\n".format(
                    name, mass, charge, ptype, *NONBONDED_CT
                )
            )
            n_fixed += 1
        else:
            out.append(line)
    NONBONDED_FILE.write_text("".join(out))
    print(f"custom-ffnonbonded.itp: retyped {n_fixed}/{len(types_to_fix)} atom types to sp3 CT")


def main():
    backbone_types = set()
    for mono, (calpha_name, cbeta_name) in MONOMERS.items():
        path = LIB_DIR / f"PAMA_{mono}.oplsaa.LigParGen.ff"
        calpha_type, cbeta_type = fix_monomer(path, calpha_name, cbeta_name)
        backbone_types.add(calpha_type)
        backbone_types.add(cbeta_type)
    fix_nonbonded(backbone_types)


if __name__ == "__main__":
    main()
