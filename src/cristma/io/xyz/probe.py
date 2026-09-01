"""Lightweight content recognition for XYZ and extXYZ sources."""

from __future__ import annotations


def probe_xyz(source: str) -> float:
    """Return conservative confidence from the first complete XYZ frame."""

    lines = source.splitlines()
    if len(lines) < 2:
        return 0.0
    try:
        atom_count = int(lines[0].strip())
    except ValueError:
        return 0.0
    if atom_count < 0 or len(lines) < atom_count + 2:
        return 0.0
    for row in lines[2 : 2 + min(atom_count, 8)]:
        tokens = row.split()
        if len(tokens) < 4:
            return 0.0
        try:
            tuple(float(value) for value in tokens[1:4])
        except ValueError:
            if "Properties=" not in lines[1]:
                return 0.0
    return 0.98 if "Properties=" in lines[1] else 0.8


__all__ = ["probe_xyz"]
