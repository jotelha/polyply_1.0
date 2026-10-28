#!/usr/bin/env python3
"""Retype the PAMA alkyl side chains beyond the first ester-adjacent carbon
to L-OPLS-AA (Siu, Pluhackova & Boeckmann, JCTC 2012), per
etc/lopls-reference/oplsaa.ff (the "#ifdef LOPLS" branches).

Scope, confirmed against the vendored reference:
  - Only two nonbonded deltas matter here: plain CH3/CH2 carbon charges
    (-0.180/-0.120 -> -0.222/-0.148) and their H charges (0.060 -> 0.074),
    with CH2-H also getting a reduced LJ epsilon (0.12552 -> 0.11000).
  - Branch-point / quaternary carbons within the tail (e.g. EH's
    2-ethylhexyl branch) have no tabulated L-OPLS delta and are left on
    standard OPLS-AA, same reasoning as the backbone quaternary carbon.
  - Only the all-carbon CT-CT-CT-CT dihedral has an L-OPLS delta that
    applies here. The ester-junction dihedrals (CT-C-OS-CT etc.) in these
    monomer files use OPLS-AA's "C" bond type, not "C_2"; the lopls
    repo's ester deltas are only defined for "C_2", so there is no L-OPLS
    value to substitute there -- left on standard OPLS-AA (a deliberate,
    documented deviation from the original plan sketch, which assumed the
    "C_2" variant applied).

Usage:
    python retype_alkyl_tail_lopls.py
"""
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
LIB_DIR = REPO_ROOT / "polyply" / "data" / "oplsaaLigParGen"
NONBONDED_FILE = REPO_ROOT / "etc" / "custom-ffnonbonded.itp"

MONOMERS = {
    "EH": ("C04", "C05"),
    "LA": ("C04", "C05"),
    "MMA": ("C07", "C08"),  # no tail; skipped
    "OC": ("C03", "C01"),
    "ST": ("C04", "C05"),
    "TD": ("C18", "C19"),
}

CH3_CHARGE, CH2_CHARGE = -0.222, -0.148
CH3_H_CHARGE, CH2_H_CHARGE = 0.074, 0.074
CH2_H_EPSILON = 0.11000  # sigma unchanged at 0.250
DIHEDRAL_CT_CT_CT_CT_LOPLS = (0.518787, -0.230192, 0.896807, -1.49134, 0.000000, 0.000000)


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


def find_first_alkoxy_c(idx_to_name, adjacency, calpha):
    heavy = [n for n in adjacency[calpha] if not is_h(idx_to_name, n)]
    carbonyl_c = next(n for n in heavy if sum(1 for m in adjacency[n] if idx_to_name[m].startswith("O")) == 2)
    ester_o = next(
        n for n in adjacency[carbonyl_c]
        if idx_to_name[n].startswith("O") and len([m for m in adjacency[n] if not is_h(idx_to_name, m)]) == 2
    )
    first_alkoxy_c = next(n for n in adjacency[ester_o] if n != carbonyl_c)
    return first_alkoxy_c, ester_o


def collect_tail(idx_to_name, adjacency, first_alkoxy_c, ester_o):
    """BFS over heavy atoms reachable from first_alkoxy_c without going
    through ester_o, plus their attached hydrogens."""
    visited = {first_alkoxy_c}
    stack = [first_alkoxy_c]
    while stack:
        cur = stack.pop()
        for nb in adjacency[cur]:
            if nb == ester_o or nb in visited:
                continue
            visited.add(nb)
            if not is_h(idx_to_name, nb):
                stack.append(nb)
    return visited


def classify_and_charge(idx_to_name, adjacency, tail, first_alkoxy_c):
    """Return {idx: (new_charge, new_epsilon_or_None)}; skips branch points
    and the first_alkoxy_c itself (handled by the library-charge pass)."""
    out = {}
    skipped_branches = []
    for idx in tail:
        if idx == first_alkoxy_c or is_h(idx_to_name, idx):
            continue
        h_neighbors = [n for n in adjacency[idx] if is_h(idx_to_name, n)]
        heavy_neighbors = [n for n in adjacency[idx] if not is_h(idx_to_name, n)]
        if len(h_neighbors) == 3 and len(heavy_neighbors) == 1:
            out[idx] = (CH3_CHARGE, None)
            for h in h_neighbors:
                out[h] = (CH3_H_CHARGE, None)
        elif len(h_neighbors) == 2 and len(heavy_neighbors) == 2:
            out[idx] = (CH2_CHARGE, None)
            for h in h_neighbors:
                out[h] = (CH2_H_CHARGE, CH2_H_EPSILON)
        else:
            skipped_branches.append(idx_to_name[idx])
    return out, skipped_branches


def apply_atom_charges(path, charge_map, idx_to_type):
    lines = path.read_text().splitlines(keepends=True)
    out, in_atoms, n = [], False, 0
    for line in lines:
        s = line.strip()
        if s.startswith("[ atoms ]"):
            in_atoms = True
            out.append(line)
            continue
        if in_atoms and s.startswith("["):
            in_atoms = False
            out.append(line)
            continue
        if in_atoms and s and not s.startswith(";"):
            t = s.split()
            nr = int(t[0])
            if nr in charge_map:
                charge, _ = charge_map[nr]
                t[6] = "{:.8f}".format(charge)
                out.append(
                    "{:>6s}{:>11s}{:>7s}{:>7s}{:>7s}{:>7s}{:>13s}{:>11s}  ; L-OPLS\n".format(*t)
                )
                n += 1
                continue
        out.append(line)
    path.write_text("".join(out))
    return n


def apply_tail_dihedrals(path, idx_to_name, tail):
    lines = path.read_text().splitlines(keepends=True)
    out, n = [], 0
    atoms_done = False
    for line in lines:
        s = line.strip()
        if not atoms_done:
            out.append(line)
            if s.startswith("[ bonds ]"):
                atoms_done = True
            continue
        t = s.split()
        is_num = bool(t) and not s.startswith(";") and not s.startswith("[") and \
            all(__import__("re").match(r"^-?[0-9.eE+-]+$", x) for x in t[: min(4, len(t))])
        if is_num and len(t) in (10, 11) and int(t[4]) == 3:
            ai, aj, ak, al = int(t[0]), int(t[1]), int(t[2]), int(t[3])
            if all(x in tail for x in (ai, aj, ak, al)) and not any(
                idx_to_name[x].startswith("H") for x in (ai, aj, ak, al)
            ):
                c = DIHEDRAL_CT_CT_CT_CT_LOPLS
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:6d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}"
                    "  ; L-OPLS CT-CT-CT-CT\n".format(ai, aj, ak, al, 3, *c)
                )
                n += 1
                continue
        out.append(line)
    path.write_text("".join(out))
    return n


def apply_nonbonded_epsilon(types_needing_ch2h_epsilon):
    text = NONBONDED_FILE.read_text()
    out, n = [], 0
    for line in text.splitlines(keepends=True):
        t = line.split()
        if t and t[0] in types_needing_ch2h_epsilon:
            name, mass, charge, ptype, sigma = t[0], float(t[1]), t[2], t[3], float(t[4])
            out.append(
                "{:<12s}{:>10.4f}    {:>5s}    {} {:.5E} {:.5E}  ; L-OPLS CH2-H epsilon\n".format(
                    name, mass, charge, ptype, sigma, CH2_H_EPSILON
                )
            )
            n += 1
        else:
            out.append(line)
    NONBONDED_FILE.write_text("".join(out))
    return n


def main():
    ch2_h_types_to_fix = set()
    for mono, (calpha_name, _cbeta_name) in MONOMERS.items():
        if mono == "MMA":
            print(f"PAMA_MMA.oplsaa.LigParGen.ff: no tail, skipped")
            continue
        path = LIB_DIR / f"PAMA_{mono}.oplsaa.LigParGen.ff"
        lines = path.read_text().splitlines()
        idx_to_name, name_to_idx, idx_to_type, adjacency = parse_structure(lines)
        calpha = name_to_idx[calpha_name]
        first_alkoxy_c, ester_o = find_first_alkoxy_c(idx_to_name, adjacency, calpha)
        tail = collect_tail(idx_to_name, adjacency, first_alkoxy_c, ester_o)
        charge_map, skipped = classify_and_charge(idx_to_name, adjacency, tail, first_alkoxy_c)

        n_charges = apply_atom_charges(path, charge_map, idx_to_type)
        n_dih = apply_tail_dihedrals(path, idx_to_name, tail)

        for idx, (charge, eps) in charge_map.items():
            if eps is not None and idx in idx_to_type:
                ch2_h_types_to_fix.add(idx_to_type[idx])

        print(
            f"{path.name}: tail_atoms={len(tail)} charges_retyped={n_charges} "
            f"dihedrals_retyped={n_dih} branch_points_left_on_OPLS={skipped}"
        )

    n_eps = apply_nonbonded_epsilon(ch2_h_types_to_fix)
    print(f"custom-ffnonbonded.itp: retyped epsilon for {n_eps} CH2-hydrogen types")


if __name__ == "__main__":
    main()
