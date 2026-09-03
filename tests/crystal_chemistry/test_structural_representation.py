from __future__ import annotations

from cristma.chemistry import InteractionLayer
from cristma.crystal_chemistry import (
    ContactClassification,
    StructuralConnection,
    StructuralConnectionKind,
    StructuralRepresentationBuilder,
    StructuralSelectionPolicy,
    StructuralUnit,
    StructuralUnitGraph,
    StructuralUnitKind,
)
from cristma.structure import PeriodicAtomRef


def unit(unit_id: str, layer: InteractionLayer) -> StructuralUnit:
    return StructuralUnit(
        unit_id=unit_id,
        kind=StructuralUnitKind.POLYHEDRON,
        atom_refs=(PeriodicAtomRef(unit_id, (0, 0, 0)),),
        source_polyhedron_id=f"polyhedron:{unit_id}",
        interaction_layers=(layer,),
        contact_classifications=(ContactClassification.PRIMARY,),
    )


def test_structural_representation_excludes_interstitial_units_and_edges() -> None:
    structural = unit("unit:MoO4", InteractionLayer.STRUCTURAL)
    interstitial = unit("unit:CaO8", InteractionLayer.INTERSTITIAL)
    connection = StructuralConnection(
        connection_id="connection:Ca-Mo",
        first_unit_id=interstitial.unit_id,
        second_unit_id=structural.unit_id,
        lattice_translation=(0, 0, 0),
        connection_kind=StructuralConnectionKind.DIRECT_CONTACT,
        interaction_layers=(InteractionLayer.INTERSTITIAL,),
        contact_classifications=(ContactClassification.PRIMARY,),
    )
    graph = StructuralUnitGraph((interstitial, structural), (connection,))
    policy = StructuralSelectionPolicy(
        included_layers=frozenset({InteractionLayer.STRUCTURAL}),
        included_classifications=frozenset({ContactClassification.PRIMARY}),
    )

    representation = StructuralRepresentationBuilder(policy).build(graph)

    assert tuple(item.unit_id for item in representation.units) == ("unit:MoO4",)
    assert representation.connections == ()
    assert representation.excluded_unit_ids == ("unit:CaO8",)
    assert representation.excluded_connection_ids == ("connection:Ca-Mo",)


def test_selection_requires_both_layer_and_geometric_classification() -> None:
    primary = unit("unit:primary", InteractionLayer.STRUCTURAL)
    secondary = StructuralUnit(
        unit_id="unit:secondary",
        kind=StructuralUnitKind.ATOM,
        atom_refs=(PeriodicAtomRef("secondary", (0, 0, 0)),),
        interaction_layers=(InteractionLayer.STRUCTURAL,),
        contact_classifications=(ContactClassification.SECONDARY,),
    )
    graph = StructuralUnitGraph((secondary, primary), ())
    policy = StructuralSelectionPolicy(
        included_layers=frozenset({InteractionLayer.STRUCTURAL}),
        included_classifications=frozenset({ContactClassification.PRIMARY}),
    )

    representation = StructuralRepresentationBuilder(policy).build(graph)

    assert tuple(item.unit_id for item in representation.units) == ("unit:primary",)
    assert representation.excluded_unit_ids == ("unit:secondary",)
