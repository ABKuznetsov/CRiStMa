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
class SpaceGroupKey:
    """Identity of one setting in the pinned 530-entry Hall catalog."""

    hall_number: int
    hall_symbol: str
    choice: str

    def __post_init__(self) -> None:
        if isinstance(self.hall_number, bool) or not 1 <= self.hall_number <= 530:
            raise ValueError("Hall number must be between 1 and 530")
        if not self.hall_symbol.strip():
            raise ValueError("Hall symbol must not be empty")


@dataclass(frozen=True, slots=True)
class SpaceGroupRecord:
    """One catalog setting and its exact operations and Wyckoff positions."""

    key: SpaceGroupKey
    number: int
    hm_short: str
    hm_full: str
    point_group: str
    centering: str
    crystal_system: str
    operations: tuple[AffineOperation, ...]
    wyckoff_positions: tuple[WyckoffPosition, ...]

    def __post_init__(self) -> None:
        if isinstance(self.number, bool) or not 1 <= self.number <= 230:
            raise ValueError("space-group number must be between 1 and 230")
        if not self.hm_short.strip() or not self.hm_full.strip():
            raise ValueError("Hermann-Mauguin symbols must not be empty")
        if not self.operations:
            raise ValueError("space-group record must contain operations")
        if not self.centering.strip() or not self.crystal_system.strip():
            raise ValueError("centering and crystal system must not be empty")
        if any(position.space_group_key != self.key for position in self.wyckoff_positions):
            raise ValueError("Wyckoff position belongs to another space-group setting")

    def definition(self, *, provenance: SymmetryProvenance) -> SpaceGroupDefinition:
        """Build the compact symmetry definition used by crystal structures."""

        return SpaceGroupDefinition(
            operations=self.operations,
            provenance=provenance,
            number=self.number,
            hm_symbol=self.hm_full,
            hall_symbol=self.key.hall_symbol,
            setting=self.key.choice or None,
        )


__all__ = ["SpaceGroupKey", "SpaceGroupRecord"]
