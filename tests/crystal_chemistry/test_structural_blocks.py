from __future__ import annotations

import pytest

from cristma.chemistry import InteractionLayer
from cristma.crystal_chemistry import (
    ContactClassification,
    PeriodicComponent,
    PeriodicConnectivityResult,
    StructuralBlockClassification,
    StructuralBlockFinder,
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


def representation() -> StructuralRepresentation:
    unit = StructuralUnit(
        unit_id="unit:A",
        kind=StructuralUnitKind.ATOM,
        atom_refs=(PeriodicAtomRef("A", (1, 0, 0)),),
        interaction_layers=(InteractionLayer.STRUCTURAL,),
        contact_classifications=(ContactClassification.PRIMARY,),
    )
    return StructuralRepresentation(
        representation_id="representation:test",
        units=(unit,),
        connections=(),
        selection_policy=POLICY,
    )


def connectivity_for_rank(rank: int) -> PeriodicConnectivityResult:
    generators = (
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
    )[:rank]
    component = PeriodicComponent(
        component_id="component:unit:A",
        unit_ids=("unit:A",),
        connection_ids=(),
        image_offsets=(("unit:A", (2, 0, 0)),),
        closure_translations=generators,
        periodic_rank=rank,
        periodic_generators=generators,
    )
    return PeriodicConnectivityResult("representation:test", (component,))


@pytest.mark.parametrize(
    ("rank", "classification"),
    (
        (0, StructuralBlockClassification.FINITE_BLOCK),
        (1, StructuralBlockClassification.ONE_PERIODIC),
        (2, StructuralBlockClassification.LAYER),
        (3, StructuralBlockClassification.FRAMEWORK),
    ),
)
def test_block_classification_follows_exact_rank(
    rank: int,
    classification: StructuralBlockClassification,
) -> None:
    result = StructuralBlockFinder().find(representation(), connectivity_for_rank(rank))

    block = result.blocks[0]
    assert block.classification is classification
    assert block.periodic_rank == rank
    assert block.atom_refs == (PeriodicAtomRef("A", (3, 0, 0)),)


def test_block_finder_rejects_connectivity_from_another_representation() -> None:
    connectivity = PeriodicConnectivityResult(
        "representation:other",
        connectivity_for_rank(0).components,
    )

    with pytest.raises(ValueError, match="representation"):
        StructuralBlockFinder().find(representation(), connectivity)
