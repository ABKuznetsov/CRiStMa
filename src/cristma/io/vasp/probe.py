"""Lightweight content recognition for native VASP structure sources."""

from __future__ import annotations


def _three_floats(line: str) -> bool:
    tokens = line.split()
    if len(tokens) < 3:
        return False
    try:
        tuple(float(value) for value in tokens[:3])
    except ValueError:
        return False
    return True


def _looks_poscar(source: str) -> bool:
    lines = source.splitlines()
    if len(lines) < 7:
        return False
    try:
        scale_count = len(tuple(float(value) for value in lines[1].split()))
    except ValueError:
        return False
    if scale_count not in {1, 3} or not all(_three_floats(lines[index]) for index in range(2, 5)):
        return False
    population_index = 5
    try:
        tuple(int(value) for value in lines[population_index].split())
    except ValueError:
        population_index += 1
        if population_index >= len(lines):
            return False
        try:
            tuple(int(value) for value in lines[population_index].split())
        except ValueError:
            return False
    mode_index = population_index + 1
    if mode_index < len(lines) and lines[mode_index].strip().lower().startswith("s"):
        mode_index += 1
    return mode_index < len(lines) and lines[mode_index].strip()[:1].lower() in {"d", "c", "k"}


def probe_vasp(source: str) -> float:
    """Return conservative confidence for POSCAR and VASP trajectory outputs."""

    prefix = source[:100_000]
    stripped = prefix.lstrip()
    if stripped.startswith("<?xml") and ("<modeling" in prefix or "<calculation" in prefix):
        return 0.98
    if "vasp." in prefix[:5_000].lower() and "POSITION" in prefix and "TOTAL-FORCE" in prefix:
        return 0.98
    if "Direct configuration=" in prefix and _looks_poscar(prefix):
        return 0.98
    if _looks_poscar(prefix):
        return 0.9
    return 0.0


__all__ = ["probe_vasp"]
