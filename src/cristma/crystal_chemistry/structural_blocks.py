"""Structural blocks projected from exact quotient-graph components."""
from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from cristma.diagnostics import Diagnostic
from .periodic_connectivity import PeriodicConnectivityResult
from .representation import StructuralRepresentation

Translation = tuple[int, int, int]


class StructuralBlockClassification(StrEnum):
    FINITE_BLOCK = "finite_block"
    ONE_PERIODIC = "one_periodic"
    LAYER = "layer"
    FRAMEWORK = "framework"


CLASSIFICATION_BY_RANK = {
    0: StructuralBlockClassification.FINITE_BLOCK,
    1: StructuralBlockClassification.ONE_PERIODIC,
    2: StructuralBlockClassification.LAYER,
    3: StructuralBlockClassification.FRAMEWORK,
}


@dataclass(frozen=True, slots=True)
class StructuralBlock:
    """One maximal connected quotient component with exact periodic rank."""
    block_id: str
    representation_id: str
    unit_orbit_ids: tuple[str, ...]
    connection_orbit_ids: tuple[str, ...]
    rank: int
    periodic_generators: tuple[Translation, ...]
    classification: StructuralBlockClassification
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        if not self.block_id or not self.representation_id or not self.unit_orbit_ids:
            raise ValueError("structural block requires identities and unit orbits")
        if self.rank not in range(4) or self.classification is not CLASSIFICATION_BY_RANK[self.rank]:
            raise ValueError("structural block classification must follow periodic rank")
        if len(self.periodic_generators) != self.rank:
            raise ValueError("structural block generator count must equal periodic rank")
        if tuple(sorted(set(self.unit_orbit_ids))) != self.unit_orbit_ids:
            raise ValueError("block unit-orbit IDs must be unique and sorted")
        if tuple(sorted(set(self.connection_orbit_ids))) != self.connection_orbit_ids:
            raise ValueError("block connection-orbit IDs must be unique and sorted")

    @property
    def periodic_rank(self) -> int:
        return self.rank


@dataclass(frozen=True, slots=True)
class StructuralBlockOrbit:
    """Compatibility wrapper: a quotient block is already one symmetry orbit."""
    block_orbit_id: str
    representative_block_id: str
    blocks: tuple[StructuralBlock, ...]
    diagnostics: tuple[Diagnostic, ...] = ()

    def __post_init__(self) -> None:
        if len(self.blocks) != 1 or self.blocks[0].block_id != self.representative_block_id:
            raise ValueError("orbit-first structural block has exactly one quotient representative")


@dataclass(frozen=True, slots=True)
class StructuralBlockResult:
    representation_id: str
    blocks: tuple[StructuralBlock, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
    block_orbits: tuple[StructuralBlockOrbit, ...] = ()

    def __post_init__(self) -> None:
        ids = tuple(x.block_id for x in self.blocks)
        if tuple(sorted(ids)) != ids or len(set(ids)) != len(ids):
            raise ValueError("structural block IDs must be unique and sorted")
        if any(x.representation_id != self.representation_id for x in self.blocks):
            raise ValueError("structural block belongs to another representation")


@dataclass(frozen=True, slots=True)
class StructuralBlockFinder:
    def get_config(self) -> dict[str, object]:
        return {}

    def clone(self, **changes: object) -> "StructuralBlockFinder":
        if changes:
            raise TypeError("unknown StructuralBlockFinder configuration: " + ", ".join(sorted(changes)))
        return self

    def find(self, representation: StructuralRepresentation,
             connectivity: PeriodicConnectivityResult, **legacy_materialization_inputs: object) -> StructuralBlockResult:
        if legacy_materialization_inputs:
            raise TypeError("structural blocks no longer accept materialized structure or atomic views")
        if connectivity.representation_id != representation.representation_id:
            raise ValueError("connectivity belongs to another structural representation")
        known_units = {x.unit_orbit_id for x in representation.unit_orbits}
        known_connections = {x.connection_orbit_id for x in representation.connection_orbits}
        seen: set[str] = set()
        blocks: list[StructuralBlock] = []
        for component in connectivity.components:
            if not set(component.unit_orbit_ids) <= known_units or not set(component.connection_orbit_ids) <= known_connections:
                raise ValueError("periodic component references an unknown graph orbit")
            if seen & set(component.unit_orbit_ids):
                raise ValueError("unit orbit occurs in more than one component")
            seen.update(component.unit_orbit_ids)
            blocks.append(StructuralBlock(
                "structural-block:" + component.component_id,
                representation.representation_id,
                component.unit_orbit_ids,
                component.connection_orbit_ids,
                component.rank,
                component.periodic_generators,
                CLASSIFICATION_BY_RANK[component.rank],
                provenance=(("method", "cristma.structural_block_finder:3"),
                            ("component_id", component.component_id)),
            ))
        if seen != known_units:
            raise ValueError("periodic connectivity does not cover every selected unit orbit")
        ordered = tuple(sorted(blocks, key=lambda x: x.block_id))
        orbits = tuple(StructuralBlockOrbit("structural-block-orbit:" + x.block_id, x.block_id, (x,)) for x in ordered)
        return StructuralBlockResult(representation.representation_id, ordered,
                                     provenance=(("method", "cristma.structural_block_finder:3"),),
                                     block_orbits=orbits)


__all__ = ["StructuralBlock", "StructuralBlockClassification", "StructuralBlockFinder",
           "StructuralBlockOrbit", "StructuralBlockResult"]
