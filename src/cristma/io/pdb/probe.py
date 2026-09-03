"""Lightweight content recognition for PDB coordinate sources."""

from __future__ import annotations


def probe_pdb(source: str) -> float:
    """Return conservative confidence for fixed-column PDB records."""

    records = {line[:6].strip().upper() for line in source.splitlines()[:200]}
    if "ATOM" in records or "HETATM" in records:
        return 0.95 if records & {"HEADER", "CRYST1", "MODEL", "END"} else 0.8
    return 0.0


__all__ = ["probe_pdb"]
