"""Canonical structural units derived from resolved crystal chemistry."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from cristma.chemistry import InteractionLayer
from cristma.diagnostics import Diagnostic
from cristma.structure import PeriodicAtomRef

from .contacts import (
    ContactClassification,
    CrystalChemistryResolution,
    ResolvedContact,
)
from .polyhedra import CoordinationPolyhedron


_ZERO_TRANSLATION = (0, 0, 0)


class StructuralUnitKind(str, Enum):
    """Scientific kind of one structural unit."""

    POLYHEDRON = "polyhedron"
    ATOM = "atom"


@dataclass(frozen=True, slots=True)
class StructuralUnit:
    """Finite unit with atom membership expressed in a local periodic frame."""

    unit_id: str
    kind: StructuralUnitKind
    atom_refs: tuple[PeriodicAtomRef, ...]
    source_contact_ids: tuple[str, ...] = ()
    source_polyhedron_id: str | None = None
    provenance: tuple[tuple[str, object], ...] = ()
    interaction_layers: tuple[InteractionLayer, ...] = ()
    contact_classifications: tuple[ContactClassification, ...] = ()

    def __post_init__(self) -> None:
        if not self.unit_id:
            raise ValueError("structural unit ID must not be empty")
        if not self.atom_refs:
            raise ValueError("structural unit must contain at least one atom reference")
        if len(set(self.atom_refs)) != len(self.atom_refs):
            raise ValueError("structural unit atom references must be unique")
        if self.kind is StructuralUnitKind.POLYHEDRON and not self.source_polyhedron_id:
            raise ValueError("polyhedron unit requires its source polyhedron ID")
        if self.kind is StructuralUnitKind.ATOM and self.source_polyhedron_id is not None:
            raise ValueError("atomic unit cannot have a source polyhedron ID")
        if len(set(self.interaction_layers)) != len(self.interaction_layers):
            raise ValueError("structural unit interaction layers must be unique")
        if len(set(self.contact_classifications)) != len(self.contact_classifications):
            raise ValueError("structural unit contact classifications must be unique")


@dataclass(frozen=True, slots=True)
class StructuralUnitBuildResult:
    """Explicit result of canonical structural-unit construction."""

    units: tuple[StructuralUnit, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()


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
    ) -> StructuralUnitBuildResult:
        known_contacts = {
            contact.geometric_contact.contact_id: contact
            for contact in resolution.contacts
        }
        units: list[StructuralUnit] = []
        represented_atom_ids: set[str] = set()

        for polyhedron in sorted(polyhedra, key=lambda item: item.polyhedron_id):
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
                provenance=(("method", "cristma.structural_unit_builder:1"),),
                interaction_layers=layers,
                contact_classifications=classifications,
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

        return StructuralUnitBuildResult(
            units=tuple(units),
            provenance=(("method", "cristma.structural_unit_builder:1"),),
        )


__all__ = [
    "StructuralUnit",
    "StructuralUnitBuildResult",
    "StructuralUnitBuilder",
    "StructuralUnitKind",
]
