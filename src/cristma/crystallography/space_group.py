"""Typed, setting-aware space-group catalog records."""

from __future__ import annotations

from dataclasses import dataclass

from cristma.symmetry.affine import AffineOperation
from cristma.symmetry.orbit import (
    SpaceGroupDefinition,
    SymmetryProvenance,
)

from .wyckoff import WyckoffPosition


@dataclass(frozen=True, slots=True)
class SpaceGroupSetting:
    """One catalog setting and its exact operations and Wyckoff positions."""

    setting_id: int
    number: int
    hall_symbol: str
    choice: str
    hm_short: str
    hm_full: str
    point_group: str
    centering: str
    crystal_system: str
    symmetry_operations: tuple[AffineOperation, ...]
    wyckoff_positions: tuple[WyckoffPosition, ...]

    def __post_init__(self) -> None:
        if isinstance(self.setting_id, bool) or not 1 <= self.setting_id <= 530:
            raise ValueError("Hall number must be between 1 and 530")
        if isinstance(self.number, bool) or not 1 <= self.number <= 230:
            raise ValueError("space-group number must be between 1 and 230")
        if not self.hall_symbol.strip():
            raise ValueError("Hall symbol must not be empty")
        if not self.hm_short.strip() or not self.hm_full.strip():
            raise ValueError("Hermann-Mauguin symbols must not be empty")
        if not self.symmetry_operations:
            raise ValueError("space-group record must contain operations")
        if not self.centering.strip() or not self.crystal_system.strip():
            raise ValueError("centering and crystal system must not be empty")
        if any(position.setting_id != self.setting_id for position in self.wyckoff_positions):
            raise ValueError("Wyckoff position belongs to another space-group setting")

    def definition(self, *, provenance: SymmetryProvenance) -> SpaceGroupDefinition:
        """Build the compact symmetry definition used by crystal structures."""

        return SpaceGroupDefinition(
            operations=self.symmetry_operations,
            provenance=provenance,
            number=self.number,
            hm_symbol=self.hm_full,
            hall_symbol=self.hall_symbol,
            setting=self.choice or None,
        )


__all__ = ["SpaceGroupSetting"]
