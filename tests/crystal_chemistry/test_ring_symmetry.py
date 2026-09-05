from __future__ import annotations

from fractions import Fraction

from cristma.chemistry import Composition, InteractionLayer
from cristma.core import MeasuredValue, UnitCell
from cristma.crystal_chemistry import (
    ContactClassification,
    StructuralConnection,
    StructuralConnectionKind,
    StructuralRepresentation,
    StructuralSelectionPolicy,
    StructuralUnit,
    StructuralUnitKind,
)
from cristma.crystal_chemistry._ring_symmetry import (
    build_ring_orbits,
    map_periodic_atom_ref,
    map_periodic_unit_ref,
)
from cristma.crystal_chemistry.rings import PeriodicUnitRef, StructuralRing
from cristma.structure import (
    AtomicPropertyTable,
    AtomicView,
    CrystalStructure,
    ExpandedAtom,
    PeriodicAtomRef,
    SiteComponent,
    SymmetryImageProvenance,
)
from cristma.symmetry import AffineOperation, SpaceGroupDefinition, parse_xyz_operation


POLICY = StructuralSelectionPolicy(
    frozenset({InteractionLayer.STRUCTURAL}),
    frozenset({ContactClassification.PRIMARY}),
)


def _atom(atom_id: str, source: str, x: float) -> ExpandedAtom:
    image = SymmetryImageProvenance("identity", (0, 0, 0))
    return ExpandedAtom(
        atom_id,
        "structure:test",
        source,
        (x, 0.0, 0.0),
        (10.0 * x, 0.0, 0.0),
        (SiteComponent("B", MeasuredValue(1.0, None, "1")),),
        None,
        image,
        (image,),
    )


def _fixture(identity_only: bool = False):
    cell = UnitCell.cubic(MeasuredValue(10.0, None, "10"))
    atoms = (
        _atom("A1", "site:A", 0.1),
        _atom("B1", "site:B", 0.2),
        _atom("C1", "site:C", 0.3),
        _atom("X1", "site:X", 0.4),
        _atom("A2", "site:A", 0.9),
        _atom("B2", "site:B", 0.8),
        _atom("C2", "site:C", 0.7),
        _atom("X2", "site:X", 0.6),
    )
    view = AtomicView(atoms, cell, (True, True, True), AtomicPropertyTable(len(atoms)))
    operations = (parse_xyz_operation("x,y,z", operation_id="identity"),)
    inversion = parse_xyz_operation("-x,y,z", operation_id="inversion")
    if not identity_only:
        operations += (inversion,)
    structure = CrystalStructure(
        "test",
        cell,
        (),
        id="structure:test",
        space_group=SpaceGroupDefinition(operations, "reported"),
    )
    units = tuple(
        StructuralUnit(
            unit_id,
            StructuralUnitKind.ATOM,
            (PeriodicAtomRef(atom_id, (0, 0, 0)),),
            interaction_layers=(InteractionLayer.STRUCTURAL,),
            contact_classifications=(ContactClassification.PRIMARY,),
        )
        for unit_id, atom_id in (
            ("UA1", "A1"), ("UB1", "B1"), ("UC1", "C1"),
            ("UA2", "A2"), ("UB2", "B2"), ("UC2", "C2"),
        )
    )

    def edge(name: str, first: str, second: str, connector: str) -> StructuralConnection:
        return StructuralConnection(
            name,
            first,
            second,
            (0, 0, 0),
            StructuralConnectionKind.SHARED_VERTEX,
            (PeriodicAtomRef(connector, (0, 0, 0)),),
            interaction_layers=(InteractionLayer.STRUCTURAL,),
            contact_classifications=(ContactClassification.PRIMARY,),
        )

    connections = (
        edge("A1-B1", "UA1", "UB1", "X1"),
        edge("B1-C1", "UB1", "UC1", "X1"),
        edge("C1-A1", "UC1", "UA1", "X1"),
        edge("A2-B2", "UA2", "UB2", "X2"),
        edge("B2-C2", "UB2", "UC2", "X2"),
        edge("C2-A2", "UC2", "UA2", "X2"),
    )
    representation = StructuralRepresentation(
        "representation:test", units, connections, POLICY
    )
    composition = Composition.from_mapping({"B": 3})
    rings = (
        StructuralRing(
            "ring:1", "block:test", representation.representation_id,
            (PeriodicUnitRef("UA1", (0, 0, 0)), PeriodicUnitRef("UB1", (0, 0, 0)), PeriodicUnitRef("UC1", (0, 0, 0))),
            ("A1-B1", "B1-C1", "C1-A1"),
            (PeriodicAtomRef("X1", (0, 0, 0)),),
            composition, (0, 0, 0),
        ),
        StructuralRing(
            "ring:2", "block:test", representation.representation_id,
            (PeriodicUnitRef("UA2", (0, 0, 0)), PeriodicUnitRef("UB2", (0, 0, 0)), PeriodicUnitRef("UC2", (0, 0, 0))),
            ("A2-B2", "B2-C2", "C2-A2"),
            (PeriodicAtomRef("X2", (0, 0, 0)),),
            composition, (0, 0, 0),
        ),
    )
    return structure, view, representation, rings, inversion


def test_symmetry_maps_periodic_atom_image_with_rotation_and_wrap_shift() -> None:
    _, view, _, _, inversion = _fixture()

    mapped = map_periodic_atom_ref(
        inversion, PeriodicAtomRef("A1", (1, 0, 0)), view
    )

    assert mapped == PeriodicAtomRef("A2", (-2, 0, 0))


def test_unit_mapping_uses_transformed_atom_membership_not_unit_name() -> None:
    _, view, representation, _, inversion = _fixture()

    mapped = map_periodic_unit_ref(
        inversion,
        PeriodicUnitRef("UA1", (0, 0, 0)),
        view,
        representation,
    )

    assert mapped == PeriodicUnitRef("UA2", (-1, 0, 0))


def test_symmetry_equivalent_rings_form_one_orbit() -> None:
    structure, view, representation, rings, _ = _fixture()

    orbits, diagnostics = build_ring_orbits(structure, view, representation, rings)

    assert diagnostics == ()
    assert len(orbits) == 1
    assert orbits[0].multiplicity == 2


def test_topology_equal_but_not_symmetry_related_rings_stay_separate() -> None:
    structure, view, representation, rings, _ = _fixture(identity_only=True)

    orbits, diagnostics = build_ring_orbits(structure, view, representation, rings)

    assert diagnostics == ()
    assert len(orbits) == 2
