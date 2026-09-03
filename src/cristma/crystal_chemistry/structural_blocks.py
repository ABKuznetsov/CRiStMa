"""Crystal-chemical blocks projected from exact periodic components."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum

from cristma.diagnostics import Diagnostic
from cristma.structure import PeriodicAtomRef

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
class StructuralBlockResult:
    """Explicit collection of blocks for one structural representation."""

    representation_id: str
    blocks: tuple[StructuralBlock, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()


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
    ) -> StructuralBlockResult:
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
        return StructuralBlockResult(
            representation_id=representation.representation_id,
            blocks=tuple(sorted(blocks, key=lambda item: item.block_id)),
            provenance=(("method", "cristma.structural_block_finder:1"),),
        )


__all__ = [
    "StructuralBlock",
    "StructuralBlockClassification",
    "StructuralBlockFinder",
    "StructuralBlockResult",
]
