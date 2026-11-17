#!/usr/bin/env python3
"""Fix within-monomer dihedrals that used OPLS-AA's generic 'C'/'O' bond-type
table rows for the ester carbonyl carbon/oxygen, instead of the 'C_2'/'O_2'
rows that opls_465/opls_466 (the carbonyl C/O's own nonbonded citations)
actually carry in stock OPLS-AA.

Three role-patterns are affected (verified by direct diff of the vendored
etc/lopls-reference/oplsaa.ff/ffbonded.itp generic-'C'/'O' rows against the
'C_2'/'O_2' rows for the same atom positions):

  C-CT-CT-HC  -> C_2-CT-CT-HC   ("aldehyde & ketone")
  CT-CT-C-O   -> CT-CT-C_2-O_2  ("aldehyde & ketone", was "propanamide")
  CT-OS-C-O   -> CT-OS-C_2-O_2  ("esters", was "benzoic acids & esters")

Every other bonded term touching the carbonyl C/O (all 6 within-monomer
angles, and the dihedrals CT-C-OS-CT, CT-CT-C-OS, HC-CT-OS-C/C-OS-CT-HC,
C-OS-CT-CT) was checked and is numerically identical whether looked up via
the generic class or via C_2/O_2 -- only these three differ. See
etc/report.md for the full comparison table.

Usage:
    python fix_carbonyl_opls_class.py
"""
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

DIH_C2_CT_CT_HC = (-0.158990, -0.476980, 0.000000, 0.635960, 0.000000, 0.000000)
DIH_CT_CT_C2_O2 = (3.106620, -3.776060, -5.137950, 5.807390, 0.000000, 0.000000)
DIH_CT_OS_C2_O2 = (21.438810, 0.000000, -21.438810, 0.000000, 0.000000, 0.000000)


def is_h(name):
    return name.startswith("H")


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


def identify_roles(idx_to_name, adjacency, calpha, cbeta):
    heavy = [n for n in adjacency[calpha] if not is_h(idx_to_name[n])]
    others = [n for n in heavy if n != cbeta]

    def oxygen_count(idx):
        return sum(1 for m in adjacency[idx] if idx_to_name[m].startswith("O"))

    carbonyl_c, alpha_methyl_c = (
        others if oxygen_count(others[0]) > oxygen_count(others[1]) else others[::-1]
    )
    carbonyl_oxygens = [n for n in adjacency[carbonyl_c] if idx_to_name[n].startswith("O")]
    ester_o, carbonyl_o = (
        carbonyl_oxygens
        if len([m for m in adjacency[carbonyl_oxygens[0]] if not is_h(idx_to_name[m])]) == 2
        else carbonyl_oxygens[::-1]
    )
    return carbonyl_c, carbonyl_o, ester_o


def classify(idx_to_name, carbonyl_c, carbonyl_o, ester_o):
    bt = {}
    for idx, name in idx_to_name.items():
        if idx == carbonyl_c:
            bt[idx] = "C"
        elif idx == carbonyl_o:
            bt[idx] = "O"
        elif idx == ester_o:
            bt[idx] = "OS"
        elif is_h(name):
            bt[idx] = "HC"
        else:
            bt[idx] = "CT"
    return bt


def fix_monomer(path, calpha_name, cbeta_name):
    lines = path.read_text().splitlines(keepends=True)
    idx_to_name, name_to_idx, adjacency = parse_structure([l.rstrip("\n") for l in lines])
    calpha, cbeta = name_to_idx[calpha_name], name_to_idx[cbeta_name]
    carbonyl_c, carbonyl_o, ester_o = identify_roles(idx_to_name, adjacency, calpha, cbeta)
    bt = classify(idx_to_name, carbonyl_c, carbonyl_o, ester_o)

    out = []
    atoms_section_done = False
    n1 = n2 = n3 = 0
    for line in lines:
        stripped = line.strip()
        if not atoms_section_done:
            out.append(line)
            if stripped.startswith("[ bonds ]"):
                atoms_section_done = True
            continue

        tokens = stripped.split(";", 1)[0].split()
        if not stripped.startswith(";") and len(tokens) == 11 and tokens[4] == "3":
            ai, aj, ak, al = (int(tokens[i]) for i in range(4))
            role = tuple(bt[x] for x in (ai, aj, ak, al))

            if role in (("C", "CT", "CT", "HC"), ("HC", "CT", "CT", "C")):
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:6d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}"
                    "  ; corrected to C_2-CT-CT-HC (carbonyl C is opls_465/C_2 class)\n".format(
                        ai, aj, ak, al, 3, *DIH_C2_CT_CT_HC
                    )
                )
                n1 += 1
                continue

            if role in (("CT", "CT", "C", "O"), ("O", "C", "CT", "CT")):
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:6d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}"
                    "  ; corrected to CT-CT-C_2-O_2 (carbonyl C/O are opls_465/466, C_2/O_2 class)\n".format(
                        ai, aj, ak, al, 3, *DIH_CT_CT_C2_O2
                    )
                )
                n2 += 1
                continue

            if role in (("CT", "OS", "C", "O"), ("O", "C", "OS", "CT")):
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:6d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}"
                    "  ; corrected to CT-OS-C_2-O_2 (carbonyl C/O are opls_465/466, C_2/O_2 class)\n".format(
                        ai, aj, ak, al, 3, *DIH_CT_OS_C2_O2
                    )
                )
                n3 += 1
                continue

        out.append(line)

    path.write_text("".join(out))
    print(f"{path.name}: C_2-CT-CT-HC={n1} CT-CT-C_2-O_2={n2} CT-OS-C_2-O_2={n3}")


def main():
    for mono, (calpha_name, cbeta_name) in MONOMERS.items():
        fix_monomer(LIB_DIR / f"PAMA_{mono}.oplsaa.LigParGen.ff", calpha_name, cbeta_name)


if __name__ == "__main__":
    main()
