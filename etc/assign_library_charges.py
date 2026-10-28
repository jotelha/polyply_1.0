#!/usr/bin/env python3
"""Replace LigParGen CM1A charges with standard OPLS-AA library charges for
the PAMA backbone + ester region (everything up to, and including, the
first ester-adjacent alkoxy carbon). The alkyl tail beyond that boundary is
left untouched here; it is retyped to L-OPLS by retype_alkyl_tail_lopls.py.

Charge template sourced from etc/lopls-reference/oplsaa.ff/atomtypes.atp /
ffnonbonded.itp (the stock, non-LOPLS OPLS-AA entries):

    opls_139  (AA C, quaternary/neopentane-type)      ->  Calpha       0.000
    opls_136  (CT, CH2)                                ->  Cbeta       -0.120
    opls_140  (HC, generic alkyl H)                    ->  Cbeta H      0.060
    opls_135  (CT, CH3)                                ->  alpha-Me C  -0.180
    opls_140                                           ->  alpha-Me H   0.060
    opls_465  (AA C: esters)                           ->  carbonyl C   0.510
    opls_466  (AA =O: esters)                          ->  carbonyl O  -0.430
    opls_467  (AA -OR: ester)                          ->  ester O     -0.330
    opls_468  (methoxy C in esters)   [R' = CH3, MMA only]  -> first-alkoxy C 0.160
    opls_469  (methoxy Hs in esters)  [R' = CH3, MMA only]  -> its H     0.030
    opls_490  (C(H2OS) ethyl ester)   [R' = CH2-..., all others] -> first-alkoxy C 0.190
    opls_140  (generic HC; no dedicated "ethyl ester CH2 H" entry found
               in the vendored reference -- documented analogy, residual
               handled by renormalize-charges.py) -> its H       0.060

See etc/report.md for the identification trail (bonded-parameter fingerprint
match against etc/lopls-reference/oplsaa.ff/ffbonded.itp).

Usage:
    python assign_library_charges.py
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

CHARGE_CALPHA = 0.000
CHARGE_CBETA = -0.120
CHARGE_ALPHA_METHYL_C = -0.180
CHARGE_GENERIC_HC = 0.060
CHARGE_CARBONYL_C = 0.510
CHARGE_CARBONYL_O = -0.430
CHARGE_ESTER_O = -0.330
CHARGE_METHOXY_C = 0.160   # R' = CH3 (opls_468)
CHARGE_METHOXY_H = 0.030   # opls_469
CHARGE_ETHYL_ESTER_C = 0.190  # R' = CH2-... (opls_490)


def parse_structure(lines):
    idx_to_name = {}
    name_to_idx = {}
    adjacency = {}
    in_atoms, in_bonds = False, False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[ atoms ]"):
            in_atoms, in_bonds = True, False
            continue
        if stripped.startswith("[ bonds ]"):
            in_atoms, in_bonds = False, True
            continue
        if stripped.startswith("[") and not stripped.startswith("[ atoms ]") and not stripped.startswith("[ bonds ]"):
            if in_bonds:
                break
            in_atoms = False
            continue
        if not stripped or stripped.startswith(";"):
            continue
        if in_atoms:
            tokens = stripped.split()
            nr, name = int(tokens[0]), tokens[4]
            idx_to_name[nr] = name
            name_to_idx[name] = nr
            adjacency[nr] = set()
        elif in_bonds:
            tokens = stripped.split()
            if len(tokens) >= 2:
                a, b = int(tokens[0]), int(tokens[1])
                adjacency.setdefault(a, set()).add(b)
                adjacency.setdefault(b, set()).add(a)
    return idx_to_name, name_to_idx, adjacency


def is_h(idx_to_name, idx):
    return idx_to_name[idx].startswith("H")


def identify_backbone_ester_roles(idx_to_name, adjacency, calpha, cbeta):
    calpha_heavy_neighbors = [n for n in adjacency[calpha] if not is_h(idx_to_name, n)]
    assert cbeta in calpha_heavy_neighbors
    others = [n for n in calpha_heavy_neighbors if n != cbeta]
    assert len(others) == 2, f"expected 2 other heavy neighbors of Calpha, got {others}"

    def oxygen_count(idx):
        return sum(1 for n in adjacency[idx] if idx_to_name[n].startswith("O"))

    carbonyl_c, alpha_methyl_c = (others if oxygen_count(others[0]) > oxygen_count(others[1]) else others[::-1])
    assert oxygen_count(carbonyl_c) == 2, "carbonyl C should have 2 oxygen neighbors"
    alpha_methyl_hs = [n for n in adjacency[alpha_methyl_c] if is_h(idx_to_name, n)]
    assert len(alpha_methyl_hs) == 3

    cbeta_hs = [n for n in adjacency[cbeta] if is_h(idx_to_name, n)]
    assert len(cbeta_hs) == 2

    oxygens = [n for n in adjacency[carbonyl_c] if idx_to_name[n].startswith("O")]
    carbonyl_o, ester_o = None, None
    for o in oxygens:
        heavy_neighbors_of_o = [n for n in adjacency[o] if not is_h(idx_to_name, n)]
        if len(heavy_neighbors_of_o) == 1:
            carbonyl_o = o
        else:
            ester_o = o
    assert carbonyl_o is not None and ester_o is not None

    first_alkoxy_c = next(n for n in adjacency[ester_o] if n != carbonyl_c)
    first_alkoxy_hs = [n for n in adjacency[first_alkoxy_c] if is_h(idx_to_name, n)]
    assert len(first_alkoxy_hs) in (2, 3), f"unexpected H count on first alkoxy C: {first_alkoxy_hs}"

    return {
        "calpha": calpha,
        "cbeta": cbeta,
        "cbeta_h": cbeta_hs,
        "alpha_methyl_c": alpha_methyl_c,
        "alpha_methyl_h": alpha_methyl_hs,
        "carbonyl_c": carbonyl_c,
        "carbonyl_o": carbonyl_o,
        "ester_o": ester_o,
        "first_alkoxy_c": first_alkoxy_c,
        "first_alkoxy_h": first_alkoxy_hs,
    }


def charge_map(roles):
    is_methoxy = len(roles["first_alkoxy_h"]) == 3  # R' = CH3 (MMA)
    m = {
        roles["calpha"]: CHARGE_CALPHA,
        roles["cbeta"]: CHARGE_CBETA,
        roles["alpha_methyl_c"]: CHARGE_ALPHA_METHYL_C,
        roles["carbonyl_c"]: CHARGE_CARBONYL_C,
        roles["carbonyl_o"]: CHARGE_CARBONYL_O,
        roles["ester_o"]: CHARGE_ESTER_O,
        roles["first_alkoxy_c"]: CHARGE_METHOXY_C if is_methoxy else CHARGE_ETHYL_ESTER_C,
    }
    for h in roles["cbeta_h"] + roles["alpha_methyl_h"]:
        m[h] = CHARGE_GENERIC_HC
    for h in roles["first_alkoxy_h"]:
        m[h] = CHARGE_METHOXY_H if is_methoxy else CHARGE_GENERIC_HC
    return m


def apply_charges(path, charges):
    lines = path.read_text().splitlines(keepends=True)
    out = []
    in_atoms = False
    n_changed = 0
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[ atoms ]"):
            in_atoms = True
            out.append(line)
            continue
        if in_atoms and stripped.startswith("["):
            in_atoms = False
            out.append(line)
            continue
        if in_atoms and stripped and not stripped.startswith(";"):
            tokens = stripped.split()
            nr = int(tokens[0])
            if nr in charges:
                tokens[6] = "{:.8f}".format(charges[nr])
                out.append(
                    "{:>6s}{:>11s}{:>7s}{:>7s}{:>7s}{:>7s}{:>13s}{:>11s}  ; OPLS-AA library charge\n".format(
                        *tokens
                    )
                )
                n_changed += 1
                continue
        out.append(line)
    path.write_text("".join(out))
    return n_changed


def main():
    for mono, (calpha_name, cbeta_name) in MONOMERS.items():
        path = LIB_DIR / f"PAMA_{mono}.oplsaa.LigParGen.ff"
        lines = path.read_text().splitlines()
        idx_to_name, name_to_idx, adjacency = parse_structure(lines)
        calpha, cbeta = name_to_idx[calpha_name], name_to_idx[cbeta_name]
        roles = identify_backbone_ester_roles(idx_to_name, adjacency, calpha, cbeta)
        charges = charge_map(roles)
        n_changed = apply_charges(path, charges)
        print(
            f"{path.name}: retyped {n_changed} charges "
            f"(firstAlkoxyC={idx_to_name[roles['first_alkoxy_c']]}, "
            f"{'methoxy/CH3' if len(roles['first_alkoxy_h']) == 3 else 'ethyl-ester/CH2'})"
        )


if __name__ == "__main__":
    main()
