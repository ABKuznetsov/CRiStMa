from __future__ import annotations

from cristma.chemistry import InteractionLayer
from cristma.core import MeasuredValue, UnitCell
from cristma.crystal_chemistry import (
    ContactClassification,
    StructuralBlock,
    StructuralBlockClassification,
    StructuralBlockResult,
    StructuralConnection,
    StructuralConnectionKind,
    StructuralRepresentation,
    StructuralSelectionPolicy,
    StructuralUnit,
    StructuralUnitKind,
)
from cristma.crystal_chemistry.ring_finder import RingFinder
from cristma.crystal_chemistry.rings import StructuralRingScope
from cristma.structure import (
    AtomicPropertyTable,
    AtomicView,
    CrystalStructure,
    ExpandedAtom,
    PeriodicAtomRef,
    SiteComponent,
    SymmetryImageProvenance,
)


POLICY = StructuralSelectionPolicy(
    included_layers=frozenset({InteractionLayer.STRUCTURAL}),
    included_classifications=frozenset({ContactClassification.PRIMARY}),
)


def _edge(
    name: str,
    first: str,
    second: str,
    connector: str | None,
    translation: tuple[int, int, int] = (0, 0, 0),
) -> StructuralConnection:
    kind = (
        StructuralConnectionKind.SHARED_VERTEX
        if connector is not None
        else StructuralConnectionKind.DIRECT_CONTACT
    )
    return StructuralConnection(
        name,
        first,
        second,
        translation,
        kind,
        (() if connector is None else (PeriodicAtomRef(connector, (0, 0, 0)),)),
        interaction_layers=(InteractionLayer.STRUCTURAL,),
        contact_classifications=(ContactClassification.PRIMARY,),
    )


def _view(elements: dict[str, str]) -> AtomicView[ExpandedAtom]:
    cell = UnitCell.cubic(MeasuredValue(10.0, None, "10"))
    atoms = tuple(
        ExpandedAtom(
            id=atom_id,
            structure_id="structure:test",
            source_site_id=f"site:{atom_id}",
            fractional=(index / 20.0, 0.0, 0.0),
            cartesian=(index / 2.0, 0.0, 0.0),
            components=(SiteComponent(element, MeasuredValue(1.0, None, "1")),),
            displacement=None,
            representative_image=SymmetryImageProvenance("identity", (0, 0, 0)),
            equivalent_images=(SymmetryImageProvenance("identity", (0, 0, 0)),),
        )
        for index, (atom_id, element) in enumerate(elements.items())
    )
    return AtomicView(atoms, cell, (True, True, True), AtomicPropertyTable(len(atoms)))


def _analysis(
    memberships: dict[str, tuple[str, ...]],
    elements: dict[str, str],
    edges: tuple[StructuralConnection, ...],
):
    view = _view(elements)
    units = tuple(
        StructuralUnit(
            unit_id,
            StructuralUnitKind.COORDINATION,
            tuple(PeriodicAtomRef(atom_id, (0, 0, 0)) for atom_id in atom_ids),
            source_coordination_id=f"coordination:{unit_id}",
            interaction_layers=(InteractionLayer.STRUCTURAL,),
            contact_classifications=(ContactClassification.PRIMARY,),
        )
        for unit_id, atom_ids in memberships.items()
    )
    representation = StructuralRepresentation(
        "representation:test", units, edges, POLICY
    )
    block = StructuralBlock(
        "block:test",
        representation.representation_id,
        tuple(memberships),
        tuple(PeriodicAtomRef(atom_id, (0, 0, 0)) for atom_id in elements),
        tuple(edge.connection_id for edge in edges),
        0,
        (),
        StructuralBlockClassification.FINITE_BLOCK,
    )
    structure = CrystalStructure.explicit("test", view.cell, (), id="structure:test")
    return structure, view, representation, StructuralBlockResult(
        representation.representation_id, (block,)
    )


def test_square_with_diagonal_keeps_triangles_not_composite_square() -> None:
    memberships = {
        "A": ("A", "Oab", "Oda", "Oac"),
        "B": ("B", "Oab", "Obc"),
        "C": ("C", "Obc", "Ocd", "Oac"),
        "D": ("D", "Ocd", "Oda"),
    }
    elements = {name: ("O" if name.startswith("O") else "B") for name in {
        atom for atoms in memberships.values() for atom in atoms
    }}
    edges = (
        _edge("AB", "A", "B", "Oab"),
        _edge("BC", "B", "C", "Obc"),
        _edge("CD", "C", "D", "Ocd"),
        _edge("DA", "D", "A", "Oda"),
        _edge("AC", "A", "C", "Oac"),
    )

    result = RingFinder().find_instances(*_analysis(memberships, elements, edges))

    assert result.status.value == "complete"
    assert {ring.size for ring in result.rings} == {3}
    assert len(result.rings) == 2


def test_nonzero_translation_walk_is_not_a_ring() -> None:
    memberships = {name: (name,) for name in ("A", "B", "C")}
    elements = {name: "Si" for name in memberships}
    edges = (
        _edge("AB", "A", "B", "A"),
        _edge("BC", "B", "C", "B"),
        _edge("CA", "C", "A", "C", (1, 0, 0)),
    )

    result = RingFinder().find_instances(*_analysis(memberships, elements, edges))

    assert result.rings == ()


def test_ring_composition_and_connectors_use_unique_periodic_atoms() -> None:
    memberships = {
        "U1": ("B1", "O1", "O2", "O4"),
        "U2": ("B2", "O2", "O3", "O5", "O6"),
        "U3": ("B3", "O3", "O1", "O7"),
    }
    elements = {
        **{f"B{index}": "B" for index in range(1, 4)},
        **{f"O{index}": "O" for index in range(1, 8)},
    }
    edges = (
        _edge("U1-U2", "U1", "U2", "O2"),
        _edge("U2-U3", "U2", "U3", "O3"),
        _edge("U3-U1", "U3", "U1", "O1"),
    )

    result = RingFinder().find_instances(*_analysis(memberships, elements, edges))

    assert len(result.rings) == 1
    ring = result.rings[0]
    assert ring.composition.normalized_formula == "B3O7"
    assert set(ring.connector_atom_refs) == {
        PeriodicAtomRef("O1", (0, 0, 0)),
        PeriodicAtomRef("O2", (0, 0, 0)),
        PeriodicAtomRef("O3", (0, 0, 0)),
    }


def test_direct_contacts_do_not_create_structural_rings() -> None:
    memberships = {name: (name,) for name in ("A", "B", "C")}
    elements = {name: "Si" for name in memberships}
    edges = (
        _edge("AB", "A", "B", None),
        _edge("BC", "B", "C", None),
        _edge("CA", "C", "A", None),
    )

    result = RingFinder().find_instances(*_analysis(memberships, elements, edges))

    assert result.rings == ()


def test_shared_small_cycle_separates_local_ring_from_framework_circuit() -> None:
    memberships = {
        "A": ("A", "Oab", "Oca", "Ofa"),
        "B": ("B", "Oab", "Obc", "Obd"),
        "C": ("C", "Obc", "Oca"),
        "D": ("D", "Obd", "Ode"),
        "E": ("E", "Ode", "Oef"),
        "F": ("F", "Oef", "Ofa"),
    }
    elements = {
        atom: ("O" if atom.startswith("O") else "B")
        for atoms in memberships.values()
        for atom in atoms
    }
    edges = (
        _edge("AB", "A", "B", "Oab"),
        _edge("BC", "B", "C", "Obc"),
        _edge("CA", "C", "A", "Oca"),
        _edge("BD", "B", "D", "Obd"),
        _edge("DE", "D", "E", "Ode"),
        _edge("EF", "E", "F", "Oef"),
        _edge("FA", "F", "A", "Ofa"),
    )

    result = RingFinder().find_instances(*_analysis(memberships, elements, edges))

    assert {(ring.size, ring.scope) for ring in result.rings} == {
        (3, StructuralRingScope.LOCAL),
        (5, StructuralRingScope.FRAMEWORK),
    }
