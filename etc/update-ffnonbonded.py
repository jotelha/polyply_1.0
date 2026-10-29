#!/usr/bin/env python3
"""update-ffnonbonded.py

Extract the [ atomtypes ] section from a newly generated library .ff file,
normalise the format (drop the optional LigParGen atom-name column), merge
with the existing entries in custom-ffnonbonded.itp, sort by opls number,
and rewrite the file in place.

Usage:
    python update-ffnonbonded.py SOURCE.ff custom-ffnonbonded.itp

Arguments:
    SOURCE.ff              Library .ff file produced by renumber-opls-types.py
                           (still contains the [ atomtypes ] section).
    custom-ffnonbonded.itp Shared nonbonded parameter file to update in place.

LigParGen sometimes emits a 7-column atomtypes line:
    opls_NNNN  ATOMNAME  MASS  CHARGE  PTYPE  SIGMA  EPSILON
The standard 6-column form expected by this library is:
    opls_NNNN  MASS  CHARGE  PTYPE  SIGMA  EPSILON
The script detects which form is present by checking whether the second token
starts with a digit (mass) or a letter (atom name to be dropped).
"""

import re
import sys
from pathlib import Path

OPLS_RE = re.compile(r'^opls_(\d+)$')
LINE_FMT = "{name:<12s}{mass:>10.4f}    {charge:>5.3f}    A {sigma:.5E} {epsilon:.5E}\n"


def parse_atomtype_line(line):
    """Parse one atomtype data line; return (opls_number, formatted_line).

    Accepts both the 6-column (no atom-name) and 7-column (with atom-name)
    LigParGen formats.  Returns None if the line is not an opls_ data line.
    """
    tokens = line.split()
    if not tokens or not OPLS_RE.match(tokens[0]):
        return None

    name = tokens[0]
    number = int(OPLS_RE.match(name).group(1))

    # Detect format: if tokens[1] starts with a digit it is the mass (6-col);
    # if it starts with a letter it is the atom name to be dropped (7-col).
    if tokens[1][0].isalpha():
        tokens = [tokens[0]] + tokens[2:]  # drop atom-name column

    if len(tokens) < 6:
        raise ValueError(f"Unexpected atomtype line format: {line!r}")

    _name, mass, charge, ptype, sigma, epsilon = tokens[:6]
    formatted = LINE_FMT.format(
        name=name,
        mass=float(mass),
        charge=float(charge),
        sigma=float(sigma),
        epsilon=float(epsilon),
    )
    return number, formatted


def extract_atomtypes(path):
    """Return dict {opls_number: formatted_line} from the [ atomtypes ] section."""
    entries = {}
    in_section = False
    for line in Path(path).read_text().splitlines(keepends=True):
        stripped = line.strip()
        if stripped == '[ atomtypes ]':
            in_section = True
            continue
        if in_section:
            if stripped.startswith('['):
                break
            if not stripped or stripped.startswith(';'):
                continue
            result = parse_atomtype_line(stripped)
            if result is not None:
                number, formatted = result
                entries[number] = formatted
    return entries


def read_ffnonbonded(path):
    """Return (header_lines, dict {opls_number: formatted_line}) from file."""
    header = []
    entries = {}
    p = Path(path)
    if not p.exists():
        return header, entries

    in_section = False
    for line in p.read_text().splitlines(keepends=True):
        stripped = line.strip()
        if stripped == '[ atomtypes ]':
            in_section = True
            header.append(line)
            continue
        if not in_section:
            header.append(line)
            continue
        # inside [ atomtypes ]
        if stripped.startswith('['):
            # a second section starts — stop (shouldn't happen in this file)
            break
        if not stripped or stripped.startswith(';'):
            header.append(line)  # keep comment/blank after the section header
            continue
        result = parse_atomtype_line(stripped)
        if result is not None:
            number, formatted = result
            entries[number] = formatted

    return header, entries


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    source_path = Path(sys.argv[1])
    target_path = Path(sys.argv[2])

    new_entries = extract_atomtypes(source_path)
    if not new_entries:
        print(f"WARNING: no [ atomtypes ] entries found in '{source_path}'. Nothing done.")
        sys.exit(0)
    print(f"Found {len(new_entries)} atomtype(s) in '{source_path.name}': "
          f"opls_{min(new_entries)} – opls_{max(new_entries)}")

    header, existing_entries = read_ffnonbonded(target_path)

    overlap = set(new_entries) & set(existing_entries)
    if overlap:
        nums = sorted(overlap)
        print(f"WARNING: {len(overlap)} type(s) already present, will be overwritten: "
              f"opls_{nums[0]} … opls_{nums[-1]}")

    merged = {**existing_entries, **new_entries}
    print(f"Total atomtypes after merge: {len(merged)} "
          f"(opls_{min(merged)} – opls_{max(merged)})")

    with target_path.open('w') as fout:
        for line in header:
            fout.write(line)
        for number in sorted(merged):
            fout.write(merged[number])

    print(f"Updated: {target_path}")


if __name__ == '__main__':
    main()
