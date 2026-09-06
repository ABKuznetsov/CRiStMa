"""Crystal-chemical blocks projected from exact periodic components."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import StrEnum

from cristma.diagnostics import Diagnostic
from cristma.structure import (
    AtomicView,
    CrystalStructure,
    ExpandedAtom,
    PeriodicAtomRef,
)

from .periodic_connectivity import PeriodicConnectivityResult
from .representation import StructuralRepresentation


Translation = tuple[int, int, int]


class StructuralBlockClassification(StrEnum):
    """Base dimensional classification determined only by periodic rank."""

    FINITE_BLOCK = "finite_block"
    ONE_PERIODIC = "one_periodic"
    LAYER = "layer"
    FRAMEWORK = "framework"


_CLASSIFICATION_BY_RANK = {
    0: StructuralBlockClassification.FINITE_BLOCK,
    1: StructuralBlockClassification.ONE_PERIODIC,
    2: StructuralBlockClassification.LAYER,
    3: StructuralBlockClassification.FRAMEWORK,
}


@dataclass(frozen=True, slots=True)
class StructuralBlock:
    """One maximal connected component in a selected representation."""

    block_id: str
    representation_id: str
    unit_ids: tuple[str, ...]
    atom_refs: tuple[PeriodicAtomRef, ...]
    connection_ids: tuple[str, ...]
    periodic_rank: int
    periodic_generators: tuple[Translation, ...]
    classification: StructuralBlockClassification
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
    block_orbit_id: str = field(default="", kw_only=True)

    def __post_init__(self) -> None:
        if not self.block_id or not self.representation_id or not self.unit_ids:
            raise ValueError("structural block requires identities and member units")
        if self.periodic_rank not in range(4):
            raise ValueError("structural block periodic rank must lie between zero and three")
        if self.classification is not _CLASSIFICATION_BY_RANK[self.periodic_rank]:
            raise ValueError("structural block classification must follow periodic rank")
        if len(self.periodic_generators) != self.periodic_rank:
            raise ValueError("structural block generator count must equal periodic rank")
        if len(set(self.atom_refs)) != len(self.atom_refs):
            raise ValueError("structural block atom references must be unique")


@dataclass(frozen=True, slots=True)
class StructuralBlockOrbit:
    """One exact space-group orbit of structural blocks."""

    block_orbit_id: str
    representative_block_id: str
    blocks: tuple[StructuralBlock, ...]
    diagnostics: tuple[Diagnostic, ...] = ()

    def __post_init__(self) -> None:
        if not self.block_orbit_id or not self.representative_block_id or not self.blocks:
            raise ValueError("structural-block orbit requires identities and members")
        ids = tuple(item.block_id for item in self.blocks)
        if ids != tuple(sorted(set(ids))):
            raise ValueError("structural-block orbit members must be unique and sorted")
        if ids[0] != self.representative_block_id:
            raise ValueError("representative structural block must be the first member")
        if any(item.block_orbit_id != self.block_orbit_id for item in self.blocks):
            raise ValueError("structural-block orbit member has another orbit ID")
        signatures = {
            (item.representation_id, item.periodic_rank, item.classification)
            for item in self.blocks
        }
        if len(signatures) != 1:
            raise ValueError("structural-block orbit members must share context and rank")


@dataclass(frozen=True, slots=True)
class StructuralBlockResult:
    """Explicit collection of blocks for one structural representation."""

    representation_id: str
    blocks: tuple[StructuralBlock, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
    block_orbits: tuple[StructuralBlockOrbit, ...] = field(default=(), kw_only=True)

    def __post_init__(self) -> None:
        ids = tuple(item.block_id for item in self.blocks)
        if len(set(ids)) != len(ids):
            raise ValueError("structural block IDs must be unique")
        if any(item.representation_id != self.representation_id for item in self.blocks):
            raise ValueError("structural block belongs to another representation")
        if self.block_orbits:
            orbit_ids = tuple(item.block_orbit_id for item in self.block_orbits)
            if len(set(orbit_ids)) != len(orbit_ids):
                raise ValueError("structural-block orbit IDs must be unique")
            observed = tuple(
                item.block_id for orbit in self.block_orbits for item in orbit.blocks
            )
            if len(set(observed)) != len(observed) or set(observed) != set(ids):
                raise ValueError("every structural block must belong to exactly one orbit")
        elif any(item.block_orbit_id for item in self.blocks):
            raise ValueError("block orbit IDs require structural-block orbit records")


def _add_translation(first: Translation, second: Translation) -> Translation:
    return (
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    )


@dataclass(frozen=True, slots=True)
class StructuralBlockFinder:
    """Project analysed components into block records without new inference."""

    def get_config(self) -> dict[str, object]:
        return {}

    def clone(self, **changes: object) -> "StructuralBlockFinder":
        if changes:
            names = ", ".join(sorted(changes))
            raise TypeError(f"unknown StructuralBlockFinder configuration: {names}")
        return self

    def find(
        self,
        representation: StructuralRepresentation,
        connectivity: PeriodicConnectivityResult,
        *,
        structure: CrystalStructure | None = None,
        atomic_view: AtomicView[ExpandedAtom] | None = None,
    ) -> StructuralBlockResult:
        if structure is not None and atomic_view is None:
            raise ValueError("structural-block symmetry grouping requires an atomic view")
        if connectivity.representation_id != representation.representation_id:
            raise ValueError("connectivity belongs to another structural representation")

        unit_by_id = {unit.unit_id: unit for unit in representation.units}
        connection_ids = {
            connection.connection_id for connection in representation.connections
        }
        seen_unit_ids: set[str] = set()
        blocks: list[StructuralBlock] = []
        for component in connectivity.components:
            component_unit_ids = set(component.unit_ids)
            if not component_unit_ids <= set(unit_by_id):
                raise ValueError("periodic component references an unknown structural unit")
            if seen_unit_ids & component_unit_ids:
                raise ValueError("structural unit occurs in more than one periodic component")
            if not set(component.connection_ids) <= connection_ids:
                raise ValueError("periodic component references an unknown connection")
            seen_unit_ids |= component_unit_ids

            offsets = dict(component.image_offsets)
            if set(offsets) != component_unit_ids:
                raise ValueError("periodic component image offsets do not cover its units")
            atom_refs = {
                PeriodicAtomRef(
                    atom_ref.atom_id,
                    _add_translation(offsets[unit_id], atom_ref.cell_translation),
                )
                for unit_id in component.unit_ids
                for atom_ref in unit_by_id[unit_id].atom_refs
            }
            blocks.append(StructuralBlock(
                block_id=f"block:{component.component_id}",
                representation_id=representation.representation_id,
                unit_ids=component.unit_ids,
                atom_refs=tuple(sorted(
                    atom_refs,
                    key=lambda item: (item.atom_id, item.cell_translation),
                )),
                connection_ids=component.connection_ids,
                periodic_rank=component.periodic_rank,
                periodic_generators=component.periodic_generators,
                classification=_CLASSIFICATION_BY_RANK[component.periodic_rank],
                provenance=(
                    ("method", "cristma.structural_block_finder:1"),
                    ("component_id", component.component_id),
                ),
            ))

        if seen_unit_ids != set(unit_by_id):
            raise ValueError("periodic connectivity does not cover every selected unit")
        completed_blocks = tuple(sorted(blocks, key=lambda item: item.block_id))
        block_orbits: tuple[StructuralBlockOrbit, ...] = ()
        diagnostics: tuple[Diagnostic, ...] = ()
        if structure is not None:
            from .hierarchy_orbits import build_block_orbits

            completed_blocks, block_orbits, diagnostics = build_block_orbits(
                structure,
                atomic_view,
                representation,
                connectivity,
                completed_blocks,
            )
        return StructuralBlockResult(
            representation_id=representation.representation_id,
            blocks=completed_blocks,
            diagnostics=diagnostics,
            provenance=((
                "method",
                "cristma.structural_block_finder:2"
                if structure is not None
                else "cristma.structural_block_finder:1",
            ),),
            block_orbits=block_orbits,
        )


__all__ = [
    "StructuralBlock",
    "StructuralBlockClassification",
    "StructuralBlockFinder",
    "StructuralBlockOrbit",
    "StructuralBlockResult",
]
