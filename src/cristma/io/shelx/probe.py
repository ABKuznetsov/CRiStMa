"""Conservative content probe for SHELX instruction documents."""

from __future__ import annotations


def probe_shelx(source: str) -> float:
    """Return confidence based on a coherent set of SHELX instructions."""

    commands: set[str] = set()
    for raw_line in source.splitlines()[:200]:
        line = raw_line.split("!", 1)[0].strip()
        if not line:
            continue
        commands.add(line.split(maxsplit=1)[0].upper())
    if "CELL" in commands and "SFAC" in commands:
        return 0.95
    if "CELL" in commands and commands.intersection({"LATT", "SYMM", "TITL"}):
        return 0.8
    if {"TITL", "END"} <= commands:
        return 0.35
    if "TITL" in commands:
        return 0.15
    return 0.0


__all__ = ["probe_shelx"]
