"""Filtered views of one scientific structural quotient graph."""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from cristma.chemistry import InteractionLayer
from cristma.crystallography import SymmetryContext
from cristma.diagnostics import Diagnostic
from .shell_orbits import ShellRole
from .structural_graph import StructuralConnectionOrbit, StructuralUnitGraph
from .structural_units import StructuralUnitOrbit


@dataclass(frozen=True, slots=True)
class StructuralSelectionPolicy:
    included_layers: frozenset[InteractionLayer]
    included_shell_roles: frozenset[ShellRole]

    def __post_init__(self) -> None:
        if not self.included_layers:
            raise ValueError("structural selection requires at least one interaction layer")
        if not self.included_shell_roles:
            raise ValueError("structural selection requires at least one shell role")


@dataclass(frozen=True, slots=True)
class StructuralRepresentation:
    representation_id: str
    unit_orbits: tuple[StructuralUnitOrbit, ...]
    connection_orbits: tuple[StructuralConnectionOrbit, ...]
    selection_policy: StructuralSelectionPolicy
    excluded_unit_orbit_ids: tuple[str, ...] = ()
    excluded_connection_orbit_ids: tuple[str, ...] = ()
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
    _symmetry_context: SymmetryContext = field(repr=False, compare=False, kw_only=True)

    def __post_init__(self) -> None:
        if not self.representation_id or self._symmetry_context is None:
            raise ValueError("structural representation requires identity and symmetry context")
        known = {x.unit_orbit_id for x in self.unit_orbits}
        if len(known) != len(self.unit_orbits):
            raise ValueError("structural representation unit IDs must be unique")
        if any(x.first_unit_orbit_id not in known or x.second_unit_orbit_id not in known for x in self.connection_orbits):
            raise ValueError("representation connection references an excluded unit")

    @property
    def units(self) -> tuple[StructuralUnitOrbit, ...]:
        return self.unit_orbits

    @property
    def connections(self) -> tuple[StructuralConnectionOrbit, ...]:
        return self.connection_orbits


def _matches(layers, roles, policy: StructuralSelectionPolicy) -> bool:
    return bool(set(layers) & policy.included_layers) and (
        not roles or bool(set(roles) & policy.included_shell_roles)
    )


@dataclass(frozen=True, slots=True)
class StructuralRepresentationBuilder:
    policy: StructuralSelectionPolicy

    def get_config(self) -> dict[str, tuple[str, ...]]:
        return {
            "included_layers": tuple(sorted(x.value for x in self.policy.included_layers)),
            "included_shell_roles": tuple(sorted(x.value for x in self.policy.included_shell_roles)),
        }

    def clone(self, **changes: object) -> "StructuralRepresentationBuilder":
        return replace(self, **changes)

    def build(self, graph: StructuralUnitGraph) -> StructuralRepresentation:
        connections = tuple(x for x in graph.connection_orbits if _matches(x.interaction_layers, x.shell_roles, self.policy))
        selected = {x.unit_orbit_id for x in graph.unit_orbits if _matches(x.interaction_layers, x.shell_roles, self.policy)}
        for connection in connections:
            selected.update((connection.first_unit_orbit_id, connection.second_unit_orbit_id))
        units = tuple(x for x in graph.unit_orbits if x.unit_orbit_id in selected)
        connection_ids = {x.connection_orbit_id for x in connections}
        config = self.get_config()
        identifier = "representation:layers=" + ",".join(config["included_layers"]) + ";roles=" + ",".join(config["included_shell_roles"])
        return StructuralRepresentation(
            identifier, units, connections, self.policy,
            tuple(x.unit_orbit_id for x in graph.unit_orbits if x.unit_orbit_id not in selected),
            tuple(x.connection_orbit_id for x in graph.connection_orbits if x.connection_orbit_id not in connection_ids),
            graph.diagnostics,
            (("method", "cristma.structural_representation_builder:2"), ("selection", config)),
            _symmetry_context=graph._symmetry_context,
        )


__all__ = ["StructuralRepresentation", "StructuralRepresentationBuilder", "StructuralSelectionPolicy"]
