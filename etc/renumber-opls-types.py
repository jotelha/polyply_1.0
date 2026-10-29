#!/usr/bin/env python3
"""renumber-opls-types.py

Renumber OPLS atom types in a LigParGen ITP file so they do not conflict
with an existing polyply force-field library.

Usage:
    python renumber-opls-types.py INPUT.itp OUTPUT.ff [EXISTING1.ff ...]

Arguments:
    INPUT.itp     LigParGen-generated ITP whose opls_XXXX types will be
                  renumbered.
    OUTPUT.ff     Destination file (will be overwritten).
    EXISTING*.ff  Any number of already-integrated library .ff files.
                  The script finds the highest opls_N across all of them
                  and assigns new numbers starting at max + 1.

If no existing library files are supplied the new types start at 1000,
matching the convention used by this library.

The renumbering is applied in a single regex pass so there is no risk of
a double-substitution even if old and new number ranges overlap.
"""

import re
import sys
from pathlib import Path

OPLS_RE = re.compile(r'\bopls_(\d+)\b')
DEFAULT_START = 1000  # first number used in this library


def find_max_opls(files):
    """Return the highest opls_N found across *files* (0 if none found)."""
    max_n = 0
    for path in files:
        text = Path(path).read_text()
        for m in OPLS_RE.finditer(text):
            n = int(m.group(1))
            if n > max_n:
                max_n = n
    return max_n


def collect_unique_types(text):
    """Return unique opls numbers from *text* in order of first appearance."""
    seen = {}
    for m in OPLS_RE.finditer(text):
        n = int(m.group(1))
        if n not in seen:
            seen[n] = len(seen)
    return sorted(seen, key=seen.__getitem__)


def build_mapping(old_types, start):
    """Map each old opls number to a new consecutive number from *start*."""
    return {old: start + i for i, old in enumerate(old_types)}


def apply_mapping(text, mapping):
    """Replace every opls_N in *text* using *mapping* in a single pass."""
    def replace(m):
        n = int(m.group(1))
        new = mapping.get(n)
        return f'opls_{new}' if new is not None else m.group(0)
    return OPLS_RE.sub(replace, text)


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)

    input_path  = Path(sys.argv[1])
    output_path = Path(sys.argv[2])
    existing    = sys.argv[3:]

    # ── Step 1: determine starting number ──────────────────────────────────
    if existing:
        max_existing = find_max_opls(existing)
        start = max_existing + 1
        print(f"Highest opls type in existing library: opls_{max_existing}")
    else:
        start = DEFAULT_START
        print(f"No existing library files given; starting at opls_{start}")
    print(f"New types will start at: opls_{start}")

    # ── Step 2: collect unique types from the input ─────────────────────────
    input_text = input_path.read_text()
    old_types  = collect_unique_types(input_text)

    if not old_types:
        print("WARNING: no opls_N types found in input – writing unchanged.")
        output_path.write_text(input_text)
        return

    print(f"\nFound {len(old_types)} unique opls type(s) in '{input_path.name}':")
    print(f"  range: opls_{min(old_types)} – opls_{max(old_types)}")

    # ── Step 3: build mapping and report ───────────────────────────────────
    mapping = build_mapping(old_types, start)

    print(f"\nRenumbering map (opls_OLD -> opls_NEW):")
    for old, new in mapping.items():
        print(f"  opls_{old:>6d}  ->  opls_{new}")

    # ── Step 4: apply and write ─────────────────────────────────────────────
    output_text = apply_mapping(input_text, mapping)
    output_path.write_text(output_text)
    print(f"\nWritten: {output_path}")
    print(f"New types occupy opls_{start} – opls_{start + len(old_types) - 1}")


if __name__ == '__main__':
    main()
