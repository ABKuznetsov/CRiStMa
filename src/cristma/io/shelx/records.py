"""Typed core SHELX instruction records."""

from __future__ import annotations

from dataclasses import dataclass

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.symmetry.affine import AffineOperation

from .document import ShelxInstructionRecord


@dataclass(frozen=True, slots=True, kw_only=True)
class ShelxCellInstruction(ShelxInstructionRecord):
    wavelength: MeasuredValue
    cell: UnitCell


@dataclass(frozen=True, slots=True, kw_only=True)
class ShelxZerrInstruction(ShelxInstructionRecord):
    formula_units: MeasuredValue
    cell_uncertainties: tuple[MeasuredValue, ...]


@dataclass(frozen=True, slots=True, kw_only=True)
class ShelxLattInstruction(ShelxInstructionRecord):
    code: int


@dataclass(frozen=True, slots=True, kw_only=True)
class ShelxSymmInstruction(ShelxInstructionRecord):
    operation: AffineOperation


@dataclass(frozen=True, slots=True, kw_only=True)
class ShelxSfacInstruction(ShelxInstructionRecord):
    entries: tuple[str, ...]


@dataclass(frozen=True, slots=True, kw_only=True)
class ShelxFvarInstruction(ShelxInstructionRecord):
    values: tuple[MeasuredValue, ...]


@dataclass(frozen=True, slots=True, kw_only=True)
class ShelxPartInstruction(ShelxInstructionRecord):
    part: int
    occupancy_code: str | None


@dataclass(frozen=True, slots=True, kw_only=True)
class ShelxResiInstruction(ShelxInstructionRecord):
    residue_number: int | None
    residue_class: str | None


@dataclass(frozen=True, slots=True, kw_only=True)
class ShelxHklfInstruction(ShelxInstructionRecord):
    code: int


@dataclass(frozen=True, slots=True, kw_only=True)
class ShelxEndInstruction(ShelxInstructionRecord):
    pass


__all__ = [
    "ShelxCellInstruction",
    "ShelxEndInstruction",
    "ShelxFvarInstruction",
    "ShelxHklfInstruction",
    "ShelxLattInstruction",
    "ShelxPartInstruction",
    "ShelxResiInstruction",
    "ShelxSfacInstruction",
    "ShelxSymmInstruction",
    "ShelxZerrInstruction",
]
