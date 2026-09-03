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
        GrammarOperation.CENTRE_LIGAND_SHELL,
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
        GrammarOperation.CENTRE_LIGAND_SHELL,
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


def test_unrepresented_contact_endpoints_get_one_atomic_unit_each() -> None:
    contact = resolved_contact("S1", "S2", (1, 0, 0), layer=InteractionLayer.INTRA_SUBSYSTEM)
    resolution = CrystalChemistryResolution((contact,), ())

    result = StructuralUnitBuilder().build(resolution, ())

    assert {
        (unit.kind, unit.atom_refs)
        for unit in result.units
    } == {
        (StructuralUnitKind.ATOM, (PeriodicAtomRef("S1", (0, 0, 0)),)),
        (StructuralUnitKind.ATOM, (PeriodicAtomRef("S2", (0, 0, 0)),)),
    }
