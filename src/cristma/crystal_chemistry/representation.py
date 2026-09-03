"""Explicit semantic selections over a structural-unit graph."""

from __future__ import annotations

from dataclasses import dataclass, replace

from cristma.chemistry import InteractionLayer
from cristma.diagnostics import Diagnostic

from .contacts import ContactClassification
from .structural_graph import StructuralConnection, StructuralUnitGraph
from .structural_units import StructuralUnit


@dataclass(frozen=True, slots=True)
class StructuralSelectionPolicy:
    """Caller-owned criteria for one reproducible structural representation."""

    included_layers: frozenset[InteractionLayer]
    included_classifications: frozenset[ContactClassification]

    def __post_init__(self) -> None:
        if not self.included_layers:
            raise ValueError("structural selection requires at least one interaction layer")
        if not self.included_classifications:
            raise ValueError("structural selection requires at least one contact classification")


@dataclass(frozen=True, slots=True)
class StructuralRepresentation:
    """Immutable selected view of a structural-unit graph."""

    representation_id: str
    units: tuple[StructuralUnit, ...]
    connections: tuple[StructuralConnection, ...]
    selection_policy: StructuralSelectionPolicy
    excluded_unit_ids: tuple[str, ...] = ()
    excluded_connection_ids: tuple[str, ...] = ()
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        if not self.representation_id:
            raise ValueError("structural representation ID must not be empty")
        unit_ids = {item.unit_id for item in self.units}
        if len(unit_ids) != len(self.units):
            raise ValueError("structural representation unit IDs must be unique")
        if any(
            connection.first_unit_id not in unit_ids
            or connection.second_unit_id not in unit_ids
            for connection in self.connections
        ):
            raise ValueError("representation connection references an excluded unit")


def _matches(
    layers: tuple[InteractionLayer, ...],
    classifications: tuple[ContactClassification, ...],
    policy: StructuralSelectionPolicy,
) -> bool:
    return bool(
        set(layers) & policy.included_layers
        and set(classifications) & policy.included_classifications
    )


@dataclass(frozen=True, slots=True)
class StructuralRepresentationBuilder:
    """Apply an explicit selection without recalculating scientific evidence."""

    policy: StructuralSelectionPolicy

    def get_config(self) -> dict[str, tuple[str, ...]]:
        return {
            "included_layers": tuple(sorted(item.value for item in self.policy.included_layers)),
            "included_classifications": tuple(sorted(
                item.value for item in self.policy.included_classifications
            )),
        }

    def clone(self, **changes: object) -> "StructuralRepresentationBuilder":
        return replace(self, **changes)

    def build(self, graph: StructuralUnitGraph) -> StructuralRepresentation:
        selected_units = tuple(
            unit for unit in graph.units
            if _matches(
                unit.interaction_layers,
                unit.contact_classifications,
                self.policy,
            )
        )
        selected_unit_ids = {item.unit_id for item in selected_units}
        selected_connections = tuple(
            connection for connection in graph.connections
            if connection.first_unit_id in selected_unit_ids
            and connection.second_unit_id in selected_unit_ids
            and _matches(
                connection.interaction_layers,
                connection.contact_classifications,
                self.policy,
            )
        )
        selected_connection_ids = {
            item.connection_id for item in selected_connections
        }
        config = self.get_config()
        representation_id = (
            "representation:"
            f"layers={','.join(config['included_layers'])};"
            f"classifications={','.join(config['included_classifications'])}"
        )
        return StructuralRepresentation(
            representation_id=representation_id,
            units=selected_units,
            connections=selected_connections,
            selection_policy=self.policy,
            excluded_unit_ids=tuple(
                unit.unit_id for unit in graph.units
                if unit.unit_id not in selected_unit_ids
            ),
            excluded_connection_ids=tuple(
                connection.connection_id for connection in graph.connections
                if connection.connection_id not in selected_connection_ids
            ),
            provenance=(
                ("method", "cristma.structural_representation_builder:1"),
                ("selection", config),
            ),
        )


__all__ = [
    "StructuralRepresentation",
    "StructuralRepresentationBuilder",
    "StructuralSelectionPolicy",
]
