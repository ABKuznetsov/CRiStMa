"""Stable structure and symmetry-derived atom identity records."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class SourceReference:
    """Location of one scientific record in its source."""

    source_name: str | None = None
    format: str | None = None
    record_id: str | None = None
    start_offset: int | None = None
    end_offset: int | None = None


@dataclass(frozen=True, slots=True)
class StructureProvenance:
    """Minimal origin chain shared by all canonical structures."""

    source: SourceReference | None = None
    parent_structure_id: str | None = None


@dataclass(frozen=True, slots=True)
class ExpandedAtomRef:
    """A symmetry-derived atom with resolvable asymmetric-unit provenance."""

    id: str
    structure_id: str | None
    fractional: tuple[float, float, float]
    source_site_id: str
    representative_operation_id: str
    equivalent_operation_ids: tuple[str, ...]
    cell_translation: tuple[int, int, int]

    @property
    def independent_site_id(self) -> str:
        """Compatibility name for the source asymmetric-unit site."""

        return self.source_site_id


ExpandedSite = ExpandedAtomRef


__all__ = ["ExpandedAtomRef", "ExpandedSite", "SourceReference", "StructureProvenance"]
