"""Canonical structural units derived from resolved crystal chemistry."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from enum import Enum

from cristma.chemistry import GrammarOperation, InteractionLayer
from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import (
    AtomicView,
    CrystalStructure,
    ExpandedAtom,
    PeriodicAtomRef,
)

from .contacts import (
    ContactClassification,
    CrystalChemistryResolution,
    ResolutionStatus,
    ResolvedContact,
)
from .polyhedra import CoordinationPolyhedron


_ZERO_TRANSLATION = (0, 0, 0)


class StructuralUnitKind(str, Enum):
    """Scientific kind of one structural unit."""

    POLYHEDRON = "polyhedron"
    COORDINATION = "coordination"
    FINITE_GROUP = "finite_group"
    ATOM = "atom"


class StructuralUnitGeometryKind(str, Enum):
    """Dimensional geometry derived from the unit's periodic atom positions."""

    POINT = "point"
    LINEAR = "linear"
    PLANAR_POLYGON = "planar_polygon"
    POLYHEDRON = "polyhedron"


@dataclass(frozen=True, slots=True)
class StructuralUnitGeometry:
    """Scientific vertex/face geometry without any rendering representation."""

    kind: StructuralUnitGeometryKind
    affine_dimension: int
    vertex_atom_refs: tuple[PeriodicAtomRef, ...]
    faces: tuple[tuple[int, ...], ...] = ()
    center_atom_ref: PeriodicAtomRef | None = None
    diagnostics: tuple[Diagnostic, ...] = ()

    def __post_init__(self) -> None:
        if self.affine_dimension not in range(4):
            raise ValueError("structural-unit affine dimension must lie from zero to three")
        expected = {
            StructuralUnitGeometryKind.POINT: 0,
            StructuralUnitGeometryKind.LINEAR: 1,
            StructuralUnitGeometryKind.PLANAR_POLYGON: 2,
            StructuralUnitGeometryKind.POLYHEDRON: 3,
        }[self.kind]
        if self.affine_dimension != expected:
            raise ValueError("structural-unit geometry kind must match affine dimension")
        if not self.vertex_atom_refs:
            raise ValueError("structural-unit geometry requires at least one vertex")
        if len(set(self.vertex_atom_refs)) != len(self.vertex_atom_refs):
            raise ValueError("structural-unit geometry vertices must be unique")
        for face in self.faces:
            if len(face) < 3 or len(set(face)) != len(face):
                raise ValueError("structural-unit face requires three unique vertices")
            if any(index < 0 or index >= len(self.vertex_atom_refs) for index in face):
                raise ValueError("structural-unit face references an unknown vertex")
        if self.kind is StructuralUnitGeometryKind.PLANAR_POLYGON and len(self.faces) != 1:
            raise ValueError("planar structural-unit geometry requires one polygon face")
        if self.kind in {
            StructuralUnitGeometryKind.POINT,
            StructuralUnitGeometryKind.LINEAR,
        } and self.faces:
            raise ValueError("zero- and one-dimensional unit geometry has no faces")


@dataclass(frozen=True, slots=True)
class StructuralUnit:
    """Finite unit with atom membership expressed in a local periodic frame."""

    unit_id: str
    kind: StructuralUnitKind
    atom_refs: tuple[PeriodicAtomRef, ...]
    source_contact_ids: tuple[str, ...] = ()
    source_polyhedron_id: str | None = None
    source_coordination_id: str | None = None
    provenance: tuple[tuple[str, object], ...] = ()
    interaction_layers: tuple[InteractionLayer, ...] = ()
    contact_classifications: tuple[ContactClassification, ...] = ()
    unit_orbit_id: str = field(default="", kw_only=True)
    geometry: StructuralUnitGeometry | None = field(default=None, kw_only=True)

    def __post_init__(self) -> None:
        if not self.unit_id:
            raise ValueError("structural unit ID must not be empty")
        if not self.atom_refs:
            raise ValueError("structural unit must contain at least one atom reference")
        if len(set(self.atom_refs)) != len(self.atom_refs):
            raise ValueError("structural unit atom references must be unique")
        if self.kind is StructuralUnitKind.POLYHEDRON and not self.source_polyhedron_id:
            raise ValueError("polyhedron unit requires its source polyhedron ID")
        if self.kind is StructuralUnitKind.COORDINATION and not self.source_coordination_id:
            raise ValueError("coordination unit requires its source coordination ID")
        if self.kind in {StructuralUnitKind.ATOM, StructuralUnitKind.FINITE_GROUP} and (
            self.source_polyhedron_id is not None
            or self.source_coordination_id is not None
        ):
            raise ValueError("atomic and finite-group units cannot reference shell geometry")
        if len(set(self.interaction_layers)) != len(self.interaction_layers):
            raise ValueError("structural unit interaction layers must be unique")
        if len(set(self.contact_classifications)) != len(self.contact_classifications):
            raise ValueError("structural unit contact classifications must be unique")
        if self.geometry is not None and not isinstance(
            self.geometry, StructuralUnitGeometry
        ):
            raise TypeError("structural unit geometry must be StructuralUnitGeometry")


@dataclass(frozen=True, slots=True)
class StructuralUnitOrbit:
    """One exact space-group orbit of structural-unit instances."""

    unit_orbit_id: str
    representative_unit_id: str
    units: tuple[StructuralUnit, ...]
    diagnostics: tuple[Diagnostic, ...] = ()

    def __post_init__(self) -> None:
        if not self.unit_orbit_id or not self.representative_unit_id or not self.units:
            raise ValueError("structural-unit orbit requires identities and members")
        ids = tuple(item.unit_id for item in self.units)
        if ids != tuple(sorted(set(ids))):
            raise ValueError("structural-unit orbit members must be unique and sorted")
        if ids[0] != self.representative_unit_id:
            raise ValueError("representative structural unit must be the first member")
        if any(item.unit_orbit_id != self.unit_orbit_id for item in self.units):
            raise ValueError("structural-unit orbit member has another orbit ID")
        if len({item.kind for item in self.units}) != 1:
            raise ValueError("structural-unit orbit members must have one kind")


@dataclass(frozen=True, slots=True)
class StructuralUnitBuildResult:
    """Explicit result of canonical structural-unit construction."""

    units: tuple[StructuralUnit, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()
    unit_orbits: tuple[StructuralUnitOrbit, ...] = field(default=(), kw_only=True)

    def __post_init__(self) -> None:
        ids = tuple(item.unit_id for item in self.units)
        if len(set(ids)) != len(ids):
            raise ValueError("structural unit IDs must be unique")
        if self.unit_orbits:
            orbit_ids = tuple(item.unit_orbit_id for item in self.unit_orbits)
            if len(set(orbit_ids)) != len(orbit_ids):
                raise ValueError("structural-unit orbit IDs must be unique")
            observed = tuple(
                item.unit_id for orbit in self.unit_orbits for item in orbit.units
            )
            if len(set(observed)) != len(observed) or set(observed) != set(ids):
                raise ValueError("every structural unit must belong to exactly one orbit")
        elif any(item.unit_orbit_id for item in self.units):
            raise ValueError("unit orbit IDs require structural-unit orbit records")


def _negated(translation: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(-value for value in translation)  # type: ignore[return-value]


def _ligand_ref(
    center_atom_id: str,
    contact: ResolvedContact,
) -> PeriodicAtomRef:
    geometric = contact.geometric_contact
    translation = geometric.cell_translation or _ZERO_TRANSLATION
    if center_atom_id == geometric.first_atom_id:
        return PeriodicAtomRef(geometric.second_atom_id, translation)
    if center_atom_id == geometric.second_atom_id:
        return PeriodicAtomRef(geometric.first_atom_id, _negated(translation))
    raise ValueError("polyhedron contact does not contain its centre")


def _contact_semantics(
    contacts: tuple[ResolvedContact, ...],
) -> tuple[tuple[InteractionLayer, ...], tuple[ContactClassification, ...]]:
    return (
        tuple(sorted({item.interaction_layer for item in contacts}, key=lambda item: item.value)),
        tuple(sorted(
            {item.contact_classification for item in contacts},
            key=lambda item: item.value,
        )),
    )


@dataclass(frozen=True, slots=True)
class StructuralUnitBuilder:
    """Build units without repeating any chemistry or geometry calculation."""

    def get_config(self) -> dict[str, object]:
        return {}

    def clone(self, **changes: object) -> "StructuralUnitBuilder":
        if changes:
            names = ", ".join(sorted(changes))
            raise TypeError(f"unknown StructuralUnitBuilder configuration: {names}")
        return self

    def build(
        self,
        resolution: CrystalChemistryResolution,
        polyhedra: tuple[CoordinationPolyhedron, ...],
        *,
        structure: CrystalStructure | None = None,
        atomic_view: AtomicView[ExpandedAtom] | None = None,
    ) -> StructuralUnitBuildResult:
        if structure is not None and atomic_view is None:
            raise ValueError("structural-unit symmetry grouping requires an atomic view")
        known_contacts = {
            contact.geometric_contact.contact_id: contact
            for contact in resolution.contacts
        }
        units: list[StructuralUnit] = []
        represented_atom_ids: set[str] = set()
        polyhedron_by_center = {
            polyhedron.center_atom_id: polyhedron
            for polyhedron in polyhedra
            if polyhedron.status is ResolutionStatus.RESOLVED
        }

        for shell in sorted(
            (
                item for item in resolution.coordination_shells
                if item.status is ResolutionStatus.RESOLVED
            ),
            key=lambda item: item.center_atom_id,
        ):
            members = [PeriodicAtomRef(shell.center_atom_id, _ZERO_TRANSLATION)]
            members.extend(
                _ligand_ref(shell.center_atom_id, contact)
                for contact in shell.contacts
            )
            unique_members = tuple(dict.fromkeys(members))
            contact_ids = tuple(sorted({
                contact.geometric_contact.contact_id for contact in shell.contacts
            }))
            layers, classifications = _contact_semantics(shell.contacts)
            polyhedron = polyhedron_by_center.pop(shell.center_atom_id, None)
            represented_atom_ids.update(item.atom_id for item in unique_members)
            units.append(StructuralUnit(
                unit_id=(
                    f"unit:{polyhedron.polyhedron_id}"
                    if polyhedron is not None
                    else f"unit:coordination:{shell.center_atom_id}"
                ),
                kind=(
                    StructuralUnitKind.POLYHEDRON
                    if polyhedron is not None
                    else StructuralUnitKind.COORDINATION
                ),
                atom_refs=unique_members,
                source_contact_ids=contact_ids,
                source_polyhedron_id=(
                    polyhedron.polyhedron_id if polyhedron is not None else None
                ),
                source_coordination_id=f"coordination:{shell.center_atom_id}",
                provenance=(("method", "cristma.structural_unit_builder:2"),),
                interaction_layers=layers,
                contact_classifications=classifications,
            ))

        # Keep accepting independently supplied polyhedra: this is useful for
        # callers restoring serialized derived results without shell objects.
        for polyhedron in sorted(
            polyhedron_by_center.values(), key=lambda item: item.polyhedron_id
        ):
            members = [PeriodicAtomRef(polyhedron.center_atom_id, _ZERO_TRANSLATION)]
            contact_ids: list[str] = []
            for contact in polyhedron.vertex_contacts:
                contact_id = contact.geometric_contact.contact_id
                if contact_id not in known_contacts:
                    raise ValueError(
                        "polyhedron contains a contact absent from crystal-chemistry resolution"
                    )
                members.append(_ligand_ref(polyhedron.center_atom_id, contact))
                contact_ids.append(contact_id)

            unique_members = tuple(dict.fromkeys(members))
            layers, classifications = _contact_semantics(polyhedron.vertex_contacts)
            represented_atom_ids.update(item.atom_id for item in unique_members)
            units.append(StructuralUnit(
                unit_id=f"unit:{polyhedron.polyhedron_id}",
                kind=StructuralUnitKind.POLYHEDRON,
                atom_refs=unique_members,
                source_contact_ids=tuple(sorted(set(contact_ids))),
                source_polyhedron_id=polyhedron.polyhedron_id,
                source_coordination_id=f"coordination:{polyhedron.center_atom_id}",
                provenance=(("method", "cristma.structural_unit_builder:2"),),
                interaction_layers=layers,
                contact_classifications=classifications,
            ))

        for contact in sorted(
            (
                item for item in resolution.contacts
                if item.interaction_type is GrammarOperation.INTRA_SUBSYSTEM_BONDS
                and item.contact_classification is ContactClassification.PRIMARY
            ),
            key=lambda item: item.geometric_contact.contact_id,
        ):
            geometric = contact.geometric_contact
            translation = geometric.cell_translation or _ZERO_TRANSLATION
            members = (
                PeriodicAtomRef(geometric.first_atom_id, _ZERO_TRANSLATION),
                PeriodicAtomRef(geometric.second_atom_id, translation),
            )
            represented_atom_ids.update(item.atom_id for item in members)
            units.append(StructuralUnit(
                unit_id=f"unit:group:{geometric.contact_id}",
                kind=StructuralUnitKind.FINITE_GROUP,
                atom_refs=members,
                source_contact_ids=(geometric.contact_id,),
                provenance=(("method", "cristma.structural_unit_builder:2"),),
                interaction_layers=(contact.interaction_layer,),
                contact_classifications=(contact.contact_classification,),
            ))

        endpoint_ids = {
            atom_id
            for contact in resolution.contacts
            for atom_id in (
                contact.geometric_contact.first_atom_id,
                contact.geometric_contact.second_atom_id,
            )
        }
        for atom_id in sorted(endpoint_ids - represented_atom_ids):
            source_contacts = tuple(
                contact
                for contact in resolution.contacts
                if atom_id in (
                    contact.geometric_contact.first_atom_id,
                    contact.geometric_contact.second_atom_id,
                )
            )
            source_ids = tuple(sorted(
                contact.geometric_contact.contact_id for contact in source_contacts
            ))
            layers, classifications = _contact_semantics(source_contacts)
            units.append(StructuralUnit(
                unit_id=f"unit:atom:{atom_id}",
                kind=StructuralUnitKind.ATOM,
                atom_refs=(PeriodicAtomRef(atom_id, _ZERO_TRANSLATION),),
                source_contact_ids=source_ids,
                provenance=(("method", "cristma.structural_unit_builder:1"),),
                interaction_layers=layers,
                contact_classifications=classifications,
            ))

        diagnostics: list[Diagnostic] = []
        completed_units = tuple(units)
        if atomic_view is not None:
            from .hierarchy_orbits import (
                UNIT_GEOMETRY_INCOMPLETE,
                build_unit_geometry,
            )

            polyhedron_by_id = {item.polyhedron_id: item for item in polyhedra}
            rows: list[StructuralUnit] = []
            for unit in completed_units:
                geometry = build_unit_geometry(unit, atomic_view, polyhedron_by_id)
                if geometry is None:
                    diagnostics.append(Diagnostic(
                        severity=Severity.WARNING,
                        code=UNIT_GEOMETRY_INCOMPLETE,
                        message=f"Geometry is not uniquely available for {unit.unit_id}",
                    ))
                rows.append(replace(unit, geometry=geometry))
            completed_units = tuple(rows)

        unit_orbits: tuple[StructuralUnitOrbit, ...] = ()
        if structure is not None:
            from .hierarchy_orbits import build_unit_orbits

            completed_units, unit_orbits, orbit_diagnostics = build_unit_orbits(
                structure,
                atomic_view,
                completed_units,
            )
            diagnostics.extend(orbit_diagnostics)
        return StructuralUnitBuildResult(
            units=completed_units,
            diagnostics=tuple(dict.fromkeys(diagnostics)),
            provenance=((
                "method",
                "cristma.structural_unit_builder:3"
                if structure is not None
                else "cristma.structural_unit_builder:2",
            ),),
            unit_orbits=unit_orbits,
        )


__all__ = [
    "StructuralUnit",
    "StructuralUnitBuildResult",
    "StructuralUnitBuilder",
    "StructuralUnitGeometry",
    "StructuralUnitGeometryKind",
    "StructuralUnitKind",
    "StructuralUnitOrbit",
]
