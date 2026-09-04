from __future__ import annotations

import pytest

from cristma.chemistry import InteractionLayer
from cristma.crystal_chemistry import (
    ContactClassification,
    PeriodicConnectivityAnalyzer,
    StructuralConnection,
    StructuralConnectionKind,
    StructuralRepresentation,
    StructuralSelectionPolicy,
    StructuralUnit,
    StructuralUnitKind,
)
from cristma.structure import PeriodicAtomRef


POLICY = StructuralSelectionPolicy(
    included_layers=frozenset({InteractionLayer.STRUCTURAL}),
    included_classifications=frozenset({ContactClassification.PRIMARY}),
)


def atomic_unit(unit_id: str) -> StructuralUnit:
    return StructuralUnit(
        unit_id=unit_id,
        kind=StructuralUnitKind.ATOM,
        atom_refs=(PeriodicAtomRef(unit_id, (0, 0, 0)),),
        interaction_layers=(InteractionLayer.STRUCTURAL,),
        contact_classifications=(ContactClassification.PRIMARY,),
    )


def edge(
    first: str,
    second: str,
    translation: tuple[int, int, int],
    index: int,
) -> StructuralConnection:
    return StructuralConnection(
        connection_id=f"edge:{index}",
        first_unit_id=first,
        second_unit_id=second,
        lattice_translation=translation,
        connection_kind=StructuralConnectionKind.DIRECT_CONTACT,
        interaction_layers=(InteractionLayer.STRUCTURAL,),
        contact_classifications=(ContactClassification.PRIMARY,),
    )


def representation(
    unit_ids: tuple[str, ...],
    edges: tuple[StructuralConnection, ...],
) -> StructuralRepresentation:
    return StructuralRepresentation(
        representation_id="representation:test",
        units=tuple(atomic_unit(item) for item in unit_ids),
        connections=edges,
        selection_policy=POLICY,
    )


def test_cross_cell_tree_edge_is_finite() -> None:
    selected = representation(("A", "B"), (edge("A", "B", (1, 0, 0), 0),))

    result = PeriodicConnectivityAnalyzer().analyze(selected)

    component = result.components[0]
    assert component.periodic_rank == 0
    assert component.periodic_generators == ()


def test_periodic_self_edge_has_rank_one_without_reducing_its_period() -> None:
    selected = representation(("A",), (edge("A", "A", (2, 0, 0), 0),))

    result = PeriodicConnectivityAnalyzer().analyze(selected)

    component = result.components[0]
    assert component.periodic_rank == 1
    assert component.periodic_generators == ((2, 0, 0),)


def test_periodic_generators_form_exact_integer_subgroup_basis() -> None:
    selected = representation(
        ("A",),
        (
            edge("A", "A", (2, 0, 0), 0),
            edge("A", "A", (3, 0, 0), 1),
        ),
    )

    component = PeriodicConnectivityAnalyzer().analyze(selected).components[0]

    assert component.periodic_rank == 1
    assert component.periodic_generators == ((1, 0, 0),)


@pytest.mark.parametrize(
    ("translations", "expected_rank"),
    (
        (((1, 0, 0),), 1),
        (((1, 0, 0), (0, 1, 0)), 2),
        (((1, 0, 0), (0, 1, 0), (0, 0, 1)), 3),
        (((2, 0, 0), (1, 0, 0), (0, 3, 0)), 2),
    ),
)
def test_periodic_rank_is_exact(
    translations: tuple[tuple[int, int, int], ...],
    expected_rank: int,
) -> None:
    selected = representation(
        ("A",),
        tuple(edge("A", "A", translation, index) for index, translation in enumerate(translations)),
    )

    result = PeriodicConnectivityAnalyzer().analyze(selected)

    assert result.components[0].periodic_rank == expected_rank


def test_connectivity_is_invariant_to_unit_order_and_reversed_edges() -> None:
    forward = representation(
        ("A", "B"),
        (
            edge("A", "B", (0, 0, 0), 0),
            edge("A", "B", (1, 0, 0), 1),
        ),
    )
    reversed_graph = representation(
        ("B", "A"),
        (
            edge("B", "A", (-1, 0, 0), 1),
            edge("B", "A", (0, 0, 0), 0),
        ),
    )

    first = PeriodicConnectivityAnalyzer().analyze(forward).components[0]
    second = PeriodicConnectivityAnalyzer().analyze(reversed_graph).components[0]

    assert first.unit_ids == second.unit_ids == ("A", "B")
    assert first.periodic_rank == second.periodic_rank == 1
    assert first.periodic_generators == second.periodic_generators == ((1, 0, 0),)
