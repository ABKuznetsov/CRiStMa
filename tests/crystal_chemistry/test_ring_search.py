from __future__ import annotations

from cristma.chemistry import InteractionLayer
from cristma.crystal_chemistry import (
    ContactClassification,
    StructuralBlock,
    StructuralBlockClassification,
    StructuralConnection,
    StructuralConnectionKind,
    StructuralRepresentation,
    StructuralSelectionPolicy,
    StructuralUnit,
    StructuralUnitKind,
)
from cristma.crystal_chemistry._ring_search import find_shortest_return_paths
from cristma.crystal_chemistry.rings import PeriodicUnitRef, RingSearchPolicy
from cristma.structure import PeriodicAtomRef


POLICY = StructuralSelectionPolicy(
    included_layers=frozenset({InteractionLayer.STRUCTURAL}),
    included_classifications=frozenset({ContactClassification.PRIMARY}),
)


def _unit(name: str) -> StructuralUnit:
    return StructuralUnit(
        name,
        StructuralUnitKind.ATOM,
        (PeriodicAtomRef(f"atom:{name}", (0, 0, 0)),),
        interaction_layers=(InteractionLayer.STRUCTURAL,),
        contact_classifications=(ContactClassification.PRIMARY,),
    )


def _edge(
    name: str,
    first: str,
    second: str,
    translation: tuple[int, int, int] = (0, 0, 0),
) -> StructuralConnection:
    return StructuralConnection(
        name,
        first,
        second,
        translation,
        StructuralConnectionKind.SHARED_VERTEX,
        (PeriodicAtomRef(f"connector:{name}", (0, 0, 0)),),
        interaction_layers=(InteractionLayer.STRUCTURAL,),
        contact_classifications=(ContactClassification.PRIMARY,),
    )


def _fixture(
    unit_ids: tuple[str, ...],
    edges: tuple[StructuralConnection, ...],
    removed_id: str,
) -> tuple[StructuralRepresentation, StructuralBlock, StructuralConnection]:
    representation = StructuralRepresentation(
        "representation:test",
        tuple(_unit(name) for name in unit_ids),
        edges,
        POLICY,
    )
    block = StructuralBlock(
        "block:test",
        representation.representation_id,
        unit_ids,
        tuple(PeriodicAtomRef(f"atom:{name}", (0, 0, 0)) for name in unit_ids),
        tuple(edge.connection_id for edge in edges),
        0,
        (),
        StructuralBlockClassification.FINITE_BLOCK,
    )
    return representation, block, next(edge for edge in edges if edge.connection_id == removed_id)


def test_return_path_targets_the_exact_periodic_image() -> None:
    fixture = _fixture(
        ("A", "B", "C"),
        (
            _edge("removed", "A", "B", (1, 0, 0)),
            _edge("A-C", "A", "C"),
            _edge("C-B", "C", "B", (1, 0, 0)),
        ),
        "removed",
    )

    result = find_shortest_return_paths(*fixture, RingSearchPolicy())

    assert result.complete
    assert len(result.paths) == 1
    assert result.paths[0].states == (
        PeriodicUnitRef("A", (0, 0, 0)),
        PeriodicUnitRef("C", (0, 0, 0)),
        PeriodicUnitRef("B", (1, 0, 0)),
    )


def test_equal_shortest_return_paths_are_all_retained() -> None:
    fixture = _fixture(
        ("A", "B", "C", "D"),
        (
            _edge("removed", "A", "B"),
            _edge("A-C", "A", "C"),
            _edge("C-B", "C", "B"),
            _edge("A-D", "A", "D"),
            _edge("D-B", "D", "B"),
        ),
        "removed",
    )

    result = find_shortest_return_paths(*fixture, RingSearchPolicy())

    assert result.complete
    assert {path.states[1].unit_id for path in result.paths} == {"C", "D"}


def test_removed_edge_keeps_parallel_edge_with_another_translation() -> None:
    fixture = _fixture(
        ("A", "B"),
        (
            _edge("removed", "A", "B", (1, 0, 0)),
            _edge("parallel", "A", "B", (0, 0, 0)),
        ),
        "removed",
    )

    result = find_shortest_return_paths(*fixture, RingSearchPolicy(maximum_ring_size=3))

    assert result.paths == ()
    assert result.complete


def test_state_budget_returns_explicit_incomplete_search() -> None:
    fixture = _fixture(
        ("A", "B", "C", "D"),
        (
            _edge("removed", "A", "B"),
            _edge("A-C", "A", "C"),
            _edge("A-D", "A", "D"),
            _edge("C-B", "C", "B"),
            _edge("D-B", "D", "B"),
        ),
        "removed",
    )

    result = find_shortest_return_paths(
        *fixture,
        RingSearchPolicy(maximum_states_per_connection=2),
    )

    assert not result.complete
    assert result.limit_name == "maximum_states_per_connection"
