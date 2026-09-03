from __future__ import annotations

import pytest

from cristma.chemistry import (
    GrammarOperation,
    InteractionLayer,
    InteractionPriority,
)
from cristma.core.values import MeasuredValue
from cristma.crystallography import GeometricContact
from cristma.crystal_chemistry import (
    ComponentPairInterpretation,
    ContactClassification,
    ResolvedContact,
    StructuralConnectionKind,
    StructuralGraphBuilder,
    StructuralUnit,
    StructuralUnitKind,
)
from cristma.structure import PeriodicAtomRef, SiteComponent


def resolved_contact(
    first: str,
    second: str,
    translation: tuple[int, int, int],
) -> ResolvedContact:
    first_component = SiteComponent("Si", MeasuredValue(1.0, None, "1"))
    second_component = SiteComponent("O", MeasuredValue(1.0, None, "1"))
    interpretation = ComponentPairInterpretation(
        first_component.species,
        second_component.species,
        1.0,
        1.0,
        2.0,
        1.0,
        1.0,
        GrammarOperation.COVALENT_NETWORK,
        InteractionLayer.STRUCTURAL,
        InteractionPriority.PRIMARY,
        ("Si",),
        ("O",),
    )
    contact = GeometricContact(
        f"contact:{first}|{second}|{','.join(map(str, translation))}",
        first,
        second,
        translation,
        2.0,
        (2.0, 0.0, 0.0),
        f"site:{first}",
        f"site:{second}",
        "test",
    )
    return ResolvedContact(
        contact,
        GrammarOperation.COVALENT_NETWORK,
        InteractionLayer.STRUCTURAL,
        InteractionPriority.PRIMARY,
        ContactClassification.PRIMARY,
        (interpretation,),
        1.0,
        1.0,
        1.0,
        (),
        (),
    )


def polyhedron_unit(unit_id: str, shared_count: int) -> StructuralUnit:
    atom_refs = (PeriodicAtomRef(f"center:{unit_id}", (0, 0, 0)),) + tuple(
        PeriodicAtomRef(f"shared:{index}", (0, 0, 0))
        for index in range(shared_count)
    )
    return StructuralUnit(
        unit_id,
        StructuralUnitKind.POLYHEDRON,
        atom_refs,
        (),
        f"polyhedron:{unit_id}",
    )


@pytest.mark.parametrize(
    ("shared_count", "kind"),
    (
        (1, StructuralConnectionKind.SHARED_VERTEX),
        (2, StructuralConnectionKind.SHARED_EDGE),
        (3, StructuralConnectionKind.SHARED_FACE),
    ),
)
def test_shared_membership_classifies_polyhedron_connection(
    shared_count: int,
    kind: StructuralConnectionKind,
) -> None:
    units = (polyhedron_unit("A", shared_count), polyhedron_unit("B", shared_count))

    graph = StructuralGraphBuilder().build(units, ())

    assert len(graph.connections) == 1
    assert graph.connections[0].connection_kind is kind
    assert len(graph.connections[0].shared_atom_refs) == shared_count


def test_periodic_reverse_direct_contacts_collapse_to_one_connection() -> None:
    units = (
        StructuralUnit(
            "unit:A",
            StructuralUnitKind.ATOM,
            (PeriodicAtomRef("A", (0, 0, 0)),),
        ),
        StructuralUnit(
            "unit:B",
            StructuralUnitKind.ATOM,
            (PeriodicAtomRef("B", (0, 0, 0)),),
        ),
    )
    forward = resolved_contact("A", "B", (1, 0, 0))
    reverse = resolved_contact("B", "A", (-1, 0, 0))

    graph = StructuralGraphBuilder().build(units, (forward, reverse))

    assert len(graph.connections) == 1
    connection = graph.connections[0]
    assert connection.connection_kind is StructuralConnectionKind.DIRECT_CONTACT
    assert connection.first_unit_id == "unit:A"
    assert connection.second_unit_id == "unit:B"
    assert connection.lattice_translation == (1, 0, 0)
    assert set(connection.source_contact_ids) == {
        forward.geometric_contact.contact_id,
        reverse.geometric_contact.contact_id,
    }
    assert connection.interaction_layers == (InteractionLayer.STRUCTURAL,)
    assert connection.contact_classifications == (ContactClassification.PRIMARY,)
