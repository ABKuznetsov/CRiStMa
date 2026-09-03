"""Loss-preserving PDB source records."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class PdbCryst1Record:
    """Unit-cell and declared space-group fields from one CRYST1 row."""

    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    space_group: str
    z: int | None
    line_number: int


@dataclass(frozen=True, slots=True)
class PdbAtomRecord:
    """One coordinate row from ATOM or HETATM."""

    record_name: str
    serial: int
    name: str
    alternate_location: str
    residue_name: str
    chain_id: str
    residue_sequence: str
    cartesian: tuple[float, float, float]
    occupancy: float
    b_iso: float | None
    element: str
    line_number: int
    model_number: int


@dataclass(frozen=True, slots=True)
class PdbDocument:
    """Original PDB text plus the records understood by CRiStMa."""

    raw_source: str
    source_name: str | None
    records: tuple[str, ...]
    cryst1: PdbCryst1Record | None
    atoms: tuple[PdbAtomRecord, ...]

    def render_preserved(self) -> str:
        return self.raw_source


__all__ = ["PdbAtomRecord", "PdbCryst1Record", "PdbDocument"]
