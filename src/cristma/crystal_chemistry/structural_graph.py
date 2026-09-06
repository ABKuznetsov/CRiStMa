"""Exact symmetry-quotient graph of structural-unit orbits."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum

from cristma.chemistry import InteractionLayer
from cristma.crystallography import PeriodicSymmetryRelation, SymmetryContext
from cristma.crystallography.symmetry_context import _digest
from cristma.diagnostics import Diagnostic
from .contact_analysis import ContactAnalysisResult
from .polyhedron_orbits import PolyhedronOrbitBuildResult
from .shell_orbits import ShellRole
from .structural_units import StructuralUnitBuilder, StructuralUnitOrbit


class StructuralConnectionKind(StrEnum):
    SHARED_VERTEX = "shared_vertex"
    SHARED_EDGE = "shared_edge"
    SHARED_FACE = "shared_face"
    DIRECT_CONTACT = "direct_contact"


@dataclass(frozen=True, slots=True)
class StructuralConnectionOrbit:
    connection_orbit_id: str
    first_unit_orbit_id: str
    second_unit_orbit_id: str
    periodic_relation: PeriodicSymmetryRelation
    relation_type: StructuralConnectionKind
    connector_site_refs: tuple[str, ...]
    interaction_layers: tuple[InteractionLayer, ...]
    shell_roles: tuple[ShellRole, ...]
    source_resolved_contact_orbit_ids: tuple[str, ...]
    interpretation_ids: tuple[str, ...]
    multiplicity_in_reference_cell: int
    provenance: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        if not self.connection_orbit_id or not self.first_unit_orbit_id or not self.second_unit_orbit_id:
            raise ValueError("structural connection orbit identities must not be empty")
        for values in (self.connector_site_refs, self.source_resolved_contact_orbit_ids, self.interpretation_ids):
            if tuple(sorted(set(values))) != values:
                raise ValueError("structural connection references must be unique and sorted")
        if not self.source_resolved_contact_orbit_ids:
            raise ValueError("structural connection requires scientific contact-orbit evidence")
        if self.multiplicity_in_reference_cell <= 0:
            raise ValueError("structural connection multiplicity must be positive")

    @property
    def connection_id(self) -> str:
        return self.connection_orbit_id

    @property
    def first_unit_id(self) -> str:
        return self.first_unit_orbit_id

    @property
    def second_unit_id(self) -> str:
        return self.second_unit_orbit_id

    @property
    def connection_kind(self) -> StructuralConnectionKind:
        return self.relation_type


StructuralConnection = StructuralConnectionOrbit


@dataclass(frozen=True, slots=True)
class StructuralUnitGraph:
    unit_orbits: tuple[StructuralUnitOrbit, ...]
    connection_orbits: tuple[StructuralConnectionOrbit, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
    _symmetry_context: SymmetryContext = field(repr=False, compare=False, kw_only=True)

    def __post_init__(self) -> None:
        unit_ids = tuple(x.unit_orbit_id for x in self.unit_orbits)
        connection_ids = tuple(x.connection_orbit_id for x in self.connection_orbits)
        if tuple(sorted(unit_ids)) != unit_ids or len(set(unit_ids)) != len(unit_ids):
            raise ValueError("structural unit-orbit IDs must be unique and sorted")
        if tuple(sorted(connection_ids)) != connection_ids or len(set(connection_ids)) != len(connection_ids):
            raise ValueError("structural connection-orbit IDs must be unique and sorted")
        if self._symmetry_context is None:
            raise ValueError("structural graph requires its exact symmetry context")
        known = set(unit_ids)
        if any(x.first_unit_orbit_id not in known or x.second_unit_orbit_id not in known for x in self.connection_orbits):
            raise ValueError("structural connection references an unknown unit orbit")

    @property
    def units(self) -> tuple[StructuralUnitOrbit, ...]:
        return self.unit_orbits

    @property
    def connections(self) -> tuple[StructuralConnectionOrbit, ...]:
        return self.connection_orbits


def _shell_roles_by_contact(result: ContactAnalysisResult) -> dict[str, set[ShellRole]]:
    incidence_by_id = {x.incidence_orbit_id: x for x in result.contact_incidence_orbits}
    roles: dict[str, set[ShellRole]] = {}
    for shell in result.coordination_shell_orbits:
        if shell.selected is None:
            continue
        for role, incidence_ids in (
            (ShellRole.PRIMARY, shell.selected.primary_incidence_ids),
            (ShellRole.SECONDARY, shell.selected.secondary_incidence_ids),
        ):
            for incidence_id in incidence_ids:
                contact_id = incidence_by_id[incidence_id].resolved_contact_orbit_id
                roles.setdefault(contact_id, set()).add(role)
    return roles


@dataclass(frozen=True, slots=True)
class StructuralGraphBuilder:
    def get_config(self) -> dict[str, object]:
        return {}

    def clone(self, **changes: object) -> "StructuralGraphBuilder":
        if changes:
            raise TypeError("unknown StructuralGraphBuilder configuration: " + ", ".join(sorted(changes)))
        return self

    def build(self, contact_result: ContactAnalysisResult,
              polyhedra: PolyhedronOrbitBuildResult) -> StructuralUnitGraph:
        unit_result = StructuralUnitBuilder().build(contact_result, polyhedra)
        unit_by_site = {x.center_independent_site_id: x for x in unit_result.unit_orbits}
        if len(unit_by_site) != len(unit_result.unit_orbits):
            raise ValueError("each independent site must own exactly one structural-unit orbit")
        geometry_by_id = {x.geometry_orbit_id: x for x in contact_result.pair_table.contact_orbits}
        roles_by_contact = _shell_roles_by_contact(contact_result)
        connections: list[StructuralConnectionOrbit] = []
        context = contact_result._symmetry_context
        for resolved in contact_result.contact_orbits:
            geometry = geometry_by_id[resolved.geometry_orbit_id]
            first = unit_by_site[geometry.first_independent_site_id]
            second = unit_by_site[geometry.second_independent_site_id]
            relation = geometry.canonical_relation
            if second.unit_orbit_id < first.unit_orbit_id:
                first, second = second, first
                relation = relation.inverse(context)
            if first.unit_orbit_id == second.unit_orbit_id:
                identity = relation.operation_key == context.identity_operation_key and relation.lattice_translation == (0, 0, 0)
                if identity:
                    continue
            interpretations = resolved.interpretations
            layers = tuple(sorted({x.interaction_layer for x in interpretations}, key=lambda x: x.value))
            contact_roles = roles_by_contact.get(resolved.resolved_contact_orbit_id, set())
            is_shared = ShellRole.PRIMARY in contact_roles
            relation_type = StructuralConnectionKind.SHARED_VERTEX if is_shared else StructuralConnectionKind.DIRECT_CONTACT
            roles = tuple(sorted(contact_roles, key=lambda item: item.value))
            connector_refs = tuple(sorted({geometry.first_independent_site_id, geometry.second_independent_site_id})) if is_shared else ()
            payload = {
                "first": first.unit_orbit_id, "second": second.unit_orbit_id,
                "relation": (relation.operation_key, relation.lattice_translation),
                "source": resolved.resolved_contact_orbit_id,
            }
            connections.append(StructuralConnectionOrbit(
                "structural-connection-orbit:" + _digest(payload), first.unit_orbit_id, second.unit_orbit_id,
                relation, relation_type, connector_refs, layers, roles,
                (resolved.resolved_contact_orbit_id,), tuple(sorted(x.interpretation_id for x in interpretations)),
                geometry.multiplicity_in_reference_cell,
                (("method", "cristma.structural_graph_builder:2"),),
            ))
        return StructuralUnitGraph(
            unit_result.unit_orbits,
            tuple(sorted(connections, key=lambda x: x.connection_orbit_id)),
            unit_result.diagnostics,
            (("method", "cristma.structural_graph_builder:2"),),
            _symmetry_context=context,
        )


__all__ = ["StructuralConnection", "StructuralConnectionKind", "StructuralConnectionOrbit",
           "StructuralGraphBuilder", "StructuralUnitGraph"]
