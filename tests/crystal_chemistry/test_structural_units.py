from __future__ import annotations

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
    CoordinationShell,
    CoordinationPolyhedron,
    CrystalChemistryResolution,
    ResolvedContact,
    StructuralUnitBuilder,
    StructuralUnitKind,
)
from cristma.structure import PeriodicAtomRef, SiteComponent


def resolved_contact(
    first: str,
    second: str,
    translation: tuple[int, int, int],
    *,
    layer: InteractionLayer = InteractionLayer.STRUCTURAL,
    operation: GrammarOperation = GrammarOperation.CENTRE_LIGAND_SHELL,
) -> ResolvedContact:
    first_component = SiteComponent("Mo", MeasuredValue(1.0, None, "1"))
    second_component = SiteComponent("O", MeasuredValue(1.0, None, "1"))
    interpretation = ComponentPairInterpretation(
        first_component.species,
        second_component.species,
        1.0,
        1.0,
        2.20,
        1.0,
        1.0,
        operation,
        layer,
        InteractionPriority.PRIMARY,
        ("Mo",),
        ("O",),
    )
    geometric = GeometricContact(
        f"contact:{first}|{second}|{','.join(map(str, translation))}",
        first,
        second,
        translation,
        2.20,
        (2.20, 0.0, 0.0),
        f"site:{first}",
        f"site:{second}",
        "test",
    )
    return ResolvedContact(
        geometric,
        operation,
        layer,
        InteractionPriority.PRIMARY,
        ContactClassification.PRIMARY,
        (interpretation,),
        1.0,
        1.0,
        1.0,
        (),
        (),
    )


def test_polyhedron_unit_retains_periodic_atom_membership_and_sources() -> None:
    contact = resolved_contact("Mo", "O", (1, 0, 0))
    polyhedron = CoordinationPolyhedron(
        "polyhedron:Mo",
        "site:Mo",
        "Mo",
        (),
        (contact,),
        ((2.20, 0.0, 0.0),),
        (),
        1.0,
        (0.0, 0.0, 0.0),
        0.0,
    )
    resolution = CrystalChemistryResolution((contact,), ())

    result = StructuralUnitBuilder().build(resolution, (polyhedron,))

    assert len(result.units) == 1
    unit = result.units[0]
    assert unit.kind is StructuralUnitKind.POLYHEDRON
    assert unit.atom_refs == (
        PeriodicAtomRef("Mo", (0, 0, 0)),
        PeriodicAtomRef("O", (1, 0, 0)),
    )
    assert unit.source_contact_ids == (contact.geometric_contact.contact_id,)
    assert unit.source_polyhedron_id == polyhedron.polyhedron_id
    assert unit.interaction_layers == (InteractionLayer.STRUCTURAL,)
    assert unit.contact_classifications == (ContactClassification.PRIMARY,)


def test_resolved_non_3d_shell_becomes_coordination_unit() -> None:
    contacts = (
        resolved_contact("B", "O1", (0, 0, 0)),
        resolved_contact("B", "O2", (0, 0, 0)),
        resolved_contact("B", "O3", (0, 0, 0)),
    )
    shell = CoordinationShell.resolved("site:B", "B", contacts)
    resolution = CrystalChemistryResolution(contacts, (shell,))

    result = StructuralUnitBuilder().build(resolution, ())

    assert len(result.units) == 1
    unit = result.units[0]
    assert unit.kind is StructuralUnitKind.COORDINATION
    assert {item.atom_id for item in unit.atom_refs} == {"B", "O1", "O2", "O3"}
    assert unit.source_coordination_id == "coordination:B"


def test_intra_subsystem_pair_becomes_finite_group() -> None:
    contact = resolved_contact(
        "S1",
        "S2",
        (1, 0, 0),
        layer=InteractionLayer.INTRA_SUBSYSTEM,
        operation=GrammarOperation.INTRA_SUBSYSTEM_BONDS,
    )
    resolution = CrystalChemistryResolution((contact,), ())

    result = StructuralUnitBuilder().build(resolution, ())

    assert len(result.units) == 1
    unit = result.units[0]
    assert unit.kind is StructuralUnitKind.FINITE_GROUP
    assert unit.atom_refs == (
        PeriodicAtomRef("S1", (0, 0, 0)),
        PeriodicAtomRef("S2", (1, 0, 0)),
    )
    assert unit.source_contact_ids == (contact.geometric_contact.contact_id,)
    assert unit.interaction_layers == (InteractionLayer.INTRA_SUBSYSTEM,)
