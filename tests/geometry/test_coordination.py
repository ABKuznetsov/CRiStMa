import pytest

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.geometry import CoordinationAnalyzer, NeighborFinder
from cristma.structure import CrystalStructure, IndependentSite, MolecularAtom, MolecularStructure, SiteComponent
from cristma.symmetry import expand_structure


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def component(element: str, occupancy: float = 1.0) -> SiteComponent:
    return SiteComponent(element, number(occupancy))


def coordination_case(case: str):
    if case == "molecule":
        atoms = (
            MolecularAtom("A", "A", (component("C"),), (0.0, 0.0, 0.0)),
            MolecularAtom("B", "B", (component("O"),), (1.0, 0.0, 0.0)),
            MolecularAtom("C", "C", (component("O"),), (-1.0, 0.0, 0.0)),
        )
        view = MolecularStructure("molecule", atoms).atomic_view()
        return view, NeighborFinder(1.1).find(view), "A", 2

    site_a = IndependentSite(
        "site:A", "A", (component("C"),), (number(0.05), number(0), number(0))
    )
    site_b = IndependentSite(
        "site:B", "B", (component("O"),), (number(0.95), number(0), number(0))
    )
    crystal = CrystalStructure.explicit(
        "crystal", UnitCell.cubic(number(10)), (site_a, site_b)
    )
    view = expand_structure(crystal)
    return view, NeighborFinder(1.1).find(view), view.atoms[0].id, 1


@pytest.mark.parametrize("case", ["molecule", "crystal"])
def test_coordination_number_comes_from_geometric_neighbors(case: str) -> None:
    view, graph, center_id, expected = coordination_case(case)

    result = CoordinationAnalyzer().analyze(view, graph)
    environment = result.by_atom(center_id)

    assert environment.center_atom_id == center_id
    assert environment.coordination_number == expected
    assert len(environment.neighbors) == expected


def test_mixed_center_remains_one_geometric_coordination_environment() -> None:
    center = MolecularAtom(
        "M", "M", (component("Ca", 0.7), component("Sr", 0.3)), (0.0, 0.0, 0.0)
    )
    atoms = (
        center,
        MolecularAtom("O1", "O1", (component("O"),), (1.0, 0.0, 0.0)),
        MolecularAtom("O2", "O2", (component("O"),), (-1.0, 0.0, 0.0)),
    )
    view = MolecularStructure("mixed", atoms).atomic_view()

    environment = CoordinationAnalyzer().analyze(
        view, NeighborFinder(1.1).find(view)
    ).by_atom("M")

    assert environment.coordination_number == 2
    assert tuple(item.element for item in environment.center_components) == ("Ca", "Sr")
