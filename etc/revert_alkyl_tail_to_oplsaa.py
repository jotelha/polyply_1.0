#!/usr/bin/env python3
"""Assign standard OPLS-AA library charges/LJ/dihedrals to the PAMA alkyl
side chains beyond the first ester-adjacent carbon.

History: a 2026-10-28 pass (etc/retype_alkyl_tail_lopls.py) retyped this
same region to L-OPLS-AA instead. That was reverted on 2026-11-12 (this
script) back to plain OPLS-AA, so the whole repeat unit -- backbone, ester,
and tail -- uses one consistent force field with library charges. See
etc/report.md section 7 for the full rationale.

Tail detection (BFS from the first alkoxy carbon, excluding the path back
through the ester oxygen; classify by H/heavy-neighbour count) is the same
logic as etc/retype_alkyl_tail_lopls.py, so this script is independent of
whether a monomer was ever L-OPLS-retyped -- it works the same way on a
freshly-added monomer as it does on these five. Target values are the same
constants already established and cited elsewhere in this library:
  - CH3 carbon charge  -> -0.180 (opls_135), assign_library_charges.py's
    CHARGE_ALPHA_METHYL_C
  - CH2 carbon charge  -> -0.120 (opls_136), assign_library_charges.py's
    CHARGE_CBETA
  - CH3-H / CH2-H charge -> 0.060 (opls_140), assign_library_charges.py's
    CHARGE_GENERIC_HC
  - CH2-H LJ epsilon -> 0.12552 (sigma unchanged at 0.250), the stock
    OPLS-AA value (L-OPLS's delta, now reverted, was 0.11000)
  - Tail CT-CT-CT-CT dihedral -> fix_backbone_hybridization.py's
    DIHEDRAL_CT_CT_CT_HC = (0.6276, 1.8828, 0.0, -2.5104), OPLS-AA's single
    "X-CT-CT-X" wildcard torsion (identical for CT-CT-CT-CT/CT-CT-CT-HC/
    HC-CT-CT-HC in the vendored ffbonded.itp)
  - Branch-point / quaternary carbons within the tail (e.g. EH's
    2-ethylhexyl branch) are left untouched, same reasoning as the backbone
    quaternary carbon: no single generic OPLS-AA CT role fits unambiguously,
    so the original LigParGen-derived bonded terms stand.

Usage:
    python revert_alkyl_tail_to_oplsaa.py
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

CH3_CHARGE, CH2_CHARGE, H_CHARGE = -0.180, -0.120, 0.060
CH2_H_EPSILON = 0.12552  # sigma unchanged at 0.250
DIHEDRAL_CT_CT_CT_CT_OPLS = (0.627600, 1.882800, 0.000000, -2.510400, 0.000000, 0.000000)


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
    """Return {idx: charge}; skips first_alkoxy_c (handled by
    assign_library_charges.py's backbone/ester pass) and branch points
    (left on standard OPLS, same as the backbone quaternary carbon)."""
    out = {}
    skipped_branches = []
    for idx in tail:
        if idx == first_alkoxy_c or is_h(idx_to_name, idx):
            continue
        h_neighbors = [n for n in adjacency[idx] if is_h(idx_to_name, n)]
        heavy_neighbors = [n for n in adjacency[idx] if not is_h(idx_to_name, n)]
        if len(h_neighbors) == 3 and len(heavy_neighbors) == 1:
            out[idx] = CH3_CHARGE
            for h in h_neighbors:
                out[h] = H_CHARGE
        elif len(h_neighbors) == 2 and len(heavy_neighbors) == 2:
            out[idx] = CH2_CHARGE
            for h in h_neighbors:
                out[h] = H_CHARGE
        else:
            skipped_branches.append(idx_to_name[idx])
    return out, skipped_branches


def apply_atom_charges(path, charge_map):
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
                t[6] = "{:.8f}".format(charge_map[nr])
                out.append(
                    "{:>6s}{:>11s}{:>7s}{:>7s}{:>7s}{:>7s}{:>13s}{:>11s}"
                    "  ; OPLS-AA library charge\n".format(*t)
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
                c = DIHEDRAL_CT_CT_CT_CT_OPLS
                out.append(
                    "{:6d}{:6d}{:6d}{:6d}{:6d}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}{:12.6f}"
                    "  ; OPLS-AA (reverted from L-OPLS)\n".format(ai, aj, ak, al, 3, *c)
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
                "{:<12s}{:>10.4f}    {:>5s}    {} {:.5E} {:.5E}  ; OPLS-AA (reverted from L-OPLS)\n".format(
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
            print("PAMA_MMA.oplsaa.LigParGen.ff: no tail, skipped")
            continue
        path = LIB_DIR / f"PAMA_{mono}.oplsaa.LigParGen.ff"
        lines = path.read_text().splitlines()
        idx_to_name, name_to_idx, idx_to_type, adjacency = parse_structure(lines)
        calpha = name_to_idx[calpha_name]
        first_alkoxy_c, ester_o = find_first_alkoxy_c(idx_to_name, adjacency, calpha)
        tail = collect_tail(idx_to_name, adjacency, first_alkoxy_c, ester_o)
        charge_map, skipped = classify_and_charge(idx_to_name, adjacency, tail, first_alkoxy_c)

        n_charges = apply_atom_charges(path, charge_map)
        n_dih = apply_tail_dihedrals(path, idx_to_name, tail)

        for idx in charge_map:
            if not is_h(idx_to_name, idx):
                continue
            heavy = next(iter(adjacency[idx]))
            if sum(1 for m in adjacency[heavy] if is_h(idx_to_name, m)) == 2 and heavy in idx_to_type:
                ch2_h_types_to_fix.add(idx_to_type[idx])

        print(
            f"{path.name}: tail_atoms={len(tail)} charges_assigned={n_charges} "
            f"dihedrals_reverted={n_dih} branch_points_left_on_OPLS={skipped}"
        )

    n_eps = apply_nonbonded_epsilon(ch2_h_types_to_fix)
    print(f"custom-ffnonbonded.itp: reverted epsilon for {n_eps} CH2-hydrogen types")


if __name__ == "__main__":
    main()
