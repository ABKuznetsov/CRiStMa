"""Finite periodic quotient graph over canonical structural units."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from itertools import combinations

from cristma.chemistry import InteractionLayer
from cristma.diagnostics import Diagnostic
from cristma.structure import PeriodicAtomRef

from .contacts import ContactClassification, ResolvedContact
from .structural_units import StructuralUnit, StructuralUnitKind


Translation = tuple[int, int, int]
_ZERO_TRANSLATION: Translation = (0, 0, 0)


class StructuralConnectionKind(StrEnum):
    """How two structural units are connected."""

    SHARED_VERTEX = "shared_vertex"
    SHARED_EDGE = "shared_edge"
    SHARED_FACE = "shared_face"
    DIRECT_CONTACT = "direct_contact"


@dataclass(frozen=True, slots=True)
class StructuralConnection:
    """One canonical relation between quotient-graph unit nodes."""

    connection_id: str
    first_unit_id: str
    second_unit_id: str
    lattice_translation: Translation
    connection_kind: StructuralConnectionKind
    shared_atom_refs: tuple[PeriodicAtomRef, ...] = ()
    source_contact_ids: tuple[str, ...] = ()
    interaction_layers: tuple[InteractionLayer, ...] = ()
    contact_classifications: tuple[ContactClassification, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()

    def __post_init__(self) -> None:
        if not self.connection_id or not self.first_unit_id or not self.second_unit_id:
            raise ValueError("structural connection identities must not be empty")
        if len(self.lattice_translation) != 3 or any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in self.lattice_translation
        ):
            raise ValueError("connection lattice translation must contain three integers")
        if len(set(self.shared_atom_refs)) != len(self.shared_atom_refs):
            raise ValueError("shared atom references must be unique")
        if (
            self.connection_kind is StructuralConnectionKind.DIRECT_CONTACT
            and self.shared_atom_refs
        ):
            raise ValueError("direct contact cannot declare shared atom membership")
        if (
            self.connection_kind is not StructuralConnectionKind.DIRECT_CONTACT
            and not self.shared_atom_refs
        ):
            raise ValueError("shared-unit connection requires shared atom membership")


@dataclass(frozen=True, slots=True)
class StructuralUnitGraph:
    """Finite nodes and translation-labelled edges of a periodic unit graph."""

    units: tuple[StructuralUnit, ...]
    connections: tuple[StructuralConnection, ...]
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()


def _negated(value: Translation) -> Translation:
    return (-value[0], -value[1], -value[2])


def _added(first: Translation, second: Translation) -> Translation:
    return (
        first[0] + second[0],
        first[1] + second[1],
        first[2] + second[2],
    )


def _subtracted(first: Translation, second: Translation) -> Translation:
    return _added(first, _negated(second))


def _canonical_relation(
    first_unit_id: str,
    second_unit_id: str,
    translation: Translation,
) -> tuple[str, str, Translation, bool]:
    forward = (first_unit_id, second_unit_id, translation)
    reverse = (second_unit_id, first_unit_id, _negated(translation))
    if reverse < forward:
        return reverse[0], reverse[1], reverse[2], True
    return forward[0], forward[1], forward[2], False


def _connection_id(
    first_unit_id: str,
    second_unit_id: str,
    translation: Translation,
    kind: StructuralConnectionKind,
    channel: str | None = None,
) -> str:
    shift = ",".join(str(value) for value in translation)
    suffix = f":{channel}" if channel else ""
    return f"connection:{kind.value}{suffix}:{first_unit_id}|{second_unit_id}|{shift}"


def _shared_kind(count: int) -> StructuralConnectionKind:
    if count == 1:
        return StructuralConnectionKind.SHARED_VERTEX
    if count == 2:
        return StructuralConnectionKind.SHARED_EDGE
    return StructuralConnectionKind.SHARED_FACE


def _contact_metadata(
    contacts: tuple[ResolvedContact, ...],
) -> tuple[tuple[str, ...], tuple[InteractionLayer, ...], tuple[ContactClassification, ...]]:
    return (
        tuple(sorted({item.geometric_contact.contact_id for item in contacts})),
        tuple(sorted({item.interaction_layer for item in contacts}, key=lambda item: item.value)),
        tuple(sorted(
            {item.contact_classification for item in contacts},
            key=lambda item: item.value,
        )),
    )


@dataclass(frozen=True, slots=True)
class StructuralGraphBuilder:
    """Derive graph relations only from membership and supplied resolved contacts."""

    def get_config(self) -> dict[str, object]:
        return {}

    def clone(self, **changes: object) -> "StructuralGraphBuilder":
        if changes:
            names = ", ".join(sorted(changes))
            raise TypeError(f"unknown StructuralGraphBuilder configuration: {names}")
        return self

    def build(
        self,
        units: tuple[StructuralUnit, ...],
        contacts: tuple[ResolvedContact, ...],
    ) -> StructuralUnitGraph:
        ordered_units = tuple(sorted(units, key=lambda item: item.unit_id))
        if len({item.unit_id for item in ordered_units}) != len(ordered_units):
            raise ValueError("structural unit IDs must be unique")

        memberships: dict[str, list[tuple[StructuralUnit, PeriodicAtomRef]]] = {}
        for unit in ordered_units:
            for atom_ref in unit.atom_refs:
                memberships.setdefault(atom_ref.atom_id, []).append((unit, atom_ref))

        contact_by_id = {
            contact.geometric_contact.contact_id: contact for contact in contacts
        }
        shared_groups: dict[
            tuple[str, str, Translation],
            set[PeriodicAtomRef],
        ] = {}
        for atom_memberships in memberships.values():
            for (first_unit, first_ref), (second_unit, second_ref) in combinations(
                atom_memberships, 2
            ):
                translation = _subtracted(
                    first_ref.cell_translation,
                    second_ref.cell_translation,
                )
                if first_unit.unit_id == second_unit.unit_id and translation == _ZERO_TRANSLATION:
                    continue
                first_id, second_id, canonical_translation, reversed_relation = (
                    _canonical_relation(
                        first_unit.unit_id,
                        second_unit.unit_id,
                        translation,
                    )
                )
                reference = second_ref if reversed_relation else first_ref
                shared_groups.setdefault(
                    (first_id, second_id, canonical_translation), set()
                ).add(reference)

        unit_by_id = {unit.unit_id: unit for unit in ordered_units}
        connections: list[StructuralConnection] = []
        for (first_id, second_id, translation), shared_refs in sorted(shared_groups.items()):
            ordered_refs = tuple(sorted(
                shared_refs,
                key=lambda item: (item.atom_id, item.cell_translation),
            ))
            relevant_ids = (
                set(unit_by_id[first_id].source_contact_ids)
                | set(unit_by_id[second_id].source_contact_ids)
            )
            shared_atom_ids = {item.atom_id for item in ordered_refs}
            source_contacts = tuple(
                contact_by_id[contact_id]
                for contact_id in sorted(relevant_ids & set(contact_by_id))
                if shared_atom_ids.intersection((
                    contact_by_id[contact_id].geometric_contact.first_atom_id,
                    contact_by_id[contact_id].geometric_contact.second_atom_id,
                ))
            )
            source_ids, _, _ = _contact_metadata(source_contacts)
            # A shared atom describes membership overlap, not a new chemical
            # contact. Its semantic channel is therefore only what both units
            # have in common. Taking the union would let an interstitial unit
            # leak into a structural-only representation (for example CaO8
            # sharing O with an isolated MoO4 tetrahedron).
            layers = tuple(sorted(
                set(unit_by_id[first_id].interaction_layers)
                & set(unit_by_id[second_id].interaction_layers),
                key=lambda item: item.value,
            ))
            classifications = tuple(sorted(
                set(unit_by_id[first_id].contact_classifications)
                & set(unit_by_id[second_id].contact_classifications),
                key=lambda item: item.value,
            ))
            kind = _shared_kind(len(ordered_refs))
            connections.append(StructuralConnection(
                _connection_id(first_id, second_id, translation, kind),
                first_id,
                second_id,
                translation,
                kind,
                ordered_refs,
                source_ids,
                layers,
                classifications,
                (("method", "cristma.structural_graph_builder:1"),),
            ))

        direct_groups: dict[
            tuple[str, str, Translation, str, str, str, str],
            list[ResolvedContact],
        ] = {}
        for contact in contacts:
            geometric = contact.geometric_contact
            contact_id = geometric.contact_id
            contact_translation = geometric.cell_translation or _ZERO_TRANSLATION
            for first_unit, first_ref in memberships.get(geometric.first_atom_id, ()):
                for second_unit, second_ref in memberships.get(geometric.second_atom_id, ()):
                    # A contact used to define either endpoint unit is internal to
                    # that unit. Membership already expresses the inter-unit
                    # relation, so emitting it again as DIRECT_CONTACT would
                    # describe the same chemistry twice.
                    if (
                        (
                            first_unit.kind is not StructuralUnitKind.ATOM
                            and contact_id in first_unit.source_contact_ids
                        )
                        or (
                            second_unit.kind is not StructuralUnitKind.ATOM
                            and contact_id in second_unit.source_contact_ids
                        )
                    ):
                        continue
                    translation = _subtracted(
                        _added(contact_translation, first_ref.cell_translation),
                        second_ref.cell_translation,
                    )
                    if first_unit.unit_id == second_unit.unit_id and translation == _ZERO_TRANSLATION:
                        continue
                    first_id, second_id, canonical_translation, _ = _canonical_relation(
                        first_unit.unit_id,
                        second_unit.unit_id,
                        translation,
                    )
                    direct_groups.setdefault(
                        (
                            first_id,
                            second_id,
                            canonical_translation,
                            contact.interaction_type.value,
                            contact.interaction_layer.value,
                            contact.grammar_priority.value,
                            contact.contact_classification.value,
                        ),
                        [],
                    ).append(contact)

        for key, source_contacts in sorted(direct_groups.items()):
            (
                first_id,
                second_id,
                translation,
                interaction_type,
                interaction_layer,
                grammar_priority,
                contact_classification,
            ) = key
            unique_contacts = tuple({
                item.geometric_contact.contact_id: item for item in source_contacts
            }[contact_id] for contact_id in sorted({
                item.geometric_contact.contact_id for item in source_contacts
            }))
            source_ids, layers, classifications = _contact_metadata(unique_contacts)
            kind = StructuralConnectionKind.DIRECT_CONTACT
            channel = ":".join((
                interaction_type,
                interaction_layer,
                grammar_priority,
                contact_classification,
            ))
            connections.append(StructuralConnection(
                _connection_id(first_id, second_id, translation, kind, channel),
                first_id,
                second_id,
                translation,
                kind,
                (),
                source_ids,
                layers,
                classifications,
                (("method", "cristma.structural_graph_builder:1"),),
            ))

        ordered_connections = tuple(sorted(
            connections,
            key=lambda item: (
                item.first_unit_id,
                item.second_unit_id,
                item.lattice_translation,
                item.connection_kind.value,
            ),
        ))
        return StructuralUnitGraph(
            ordered_units,
            ordered_connections,
            provenance=(("method", "cristma.structural_graph_builder:1"),),
        )


__all__ = [
    "StructuralConnection",
    "StructuralConnectionKind",
    "StructuralGraphBuilder",
    "StructuralUnitGraph",
]
