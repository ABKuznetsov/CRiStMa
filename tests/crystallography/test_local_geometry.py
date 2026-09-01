from __future__ import annotations

import pytest

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.crystallography import geometric_contacts
from cristma.geometry import NeighborFinder
from cristma.structure import (
    CrystalStructure,
    IndependentSite,
    MolecularAtom,
    MolecularStructure,
    SiteComponent,
)


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def component(symbol: str = "C") -> SiteComponent:
    return SiteComponent(symbol, number(1.0))


def molecular_view(*rows: tuple[str, tuple[float, float, float]]):
    atoms = tuple(
        MolecularAtom(atom_id, atom_id, (component(),), cartesian)
        for atom_id, cartesian in rows
    )
    return MolecularStructure("finite", atoms).atomic_view()


def site(label: str, fractional: tuple[float, float, float]) -> IndependentSite:
    return IndependentSite(
        id=f"site:{label}",
        label=label,
        components=(component(),),
        fractional=tuple(number(value) for value in fractional),
    )


def test_finite_reverse_edges_become_one_contact() -> None:
    view = molecular_view(
        ("O", (0.0, 0.0, 0.0)),
        ("H1", (1.0, 0.0, 0.0)),
        ("H2", (-1.0, 0.0, 0.0)),
    )
    graph = NeighborFinder(cutoff=1.1).find(view)

    contacts = geometric_contacts(view, graph)

    assert len(contacts) == 2
    assert {
        (contact.first_atom_id, contact.second_atom_id)
        for contact in contacts
    } == {("H1", "O"), ("H2", "O")}
    assert all(contact.cell_translation is None for contact in contacts)
    assert all(contact.distance == pytest.approx(1.0) for contact in contacts)


def test_periodic_reverse_edges_become_one_canonical_contact() -> None:
    crystal = CrystalStructure.explicit(
        "boundary",
        UnitCell.cubic(number(10.0)),
        (site("A", (0.05, 0.0, 0.0)), site("B", (0.95, 0.0, 0.0))),
        id="structure:boundary",
    )
    view = crystal.atomic_view()
    graph = NeighborFinder(cutoff=1.1).find(view)

    contacts = geometric_contacts(view, graph)

    assert len(contacts) == 1
    contact = contacts[0]
    assert contact.cell_translation == (-1, 0, 0)
    assert contact.distance == pytest.approx(1.0)
    assert contact.first_source_site_id == "site:A"
    assert contact.second_source_site_id == "site:B"
    assert contact.vector_cartesian == pytest.approx((-1.0, 0.0, 0.0))


def test_periodic_self_images_are_not_folded_into_one_lattice_direction() -> None:
    crystal = CrystalStructure.explicit(
        "self",
        UnitCell.cubic(number(2.0)),
        (site("A", (0.0, 0.0, 0.0)),),
    )
    view = crystal.atomic_view()
    graph = NeighborFinder(cutoff=2.1).find(view)

    contacts = geometric_contacts(view, graph)

    assert len(contacts) == 3
    assert {contact.cell_translation for contact in contacts} == {
        (-1, 0, 0),
        (0, -1, 0),
        (0, 0, -1),
    }


def test_contact_ids_and_order_are_deterministic() -> None:
    view = molecular_view(
        ("B", (1.0, 0.0, 0.0)),
        ("A", (0.0, 0.0, 0.0)),
        ("C", (2.0, 0.0, 0.0)),
    )
    graph = NeighborFinder(cutoff=1.1).find(view)

    first = geometric_contacts(view, graph)
    second = geometric_contacts(view, graph)

    assert first == second
    assert tuple(item.contact_id for item in first) == tuple(
        sorted(item.contact_id for item in first)
    )
