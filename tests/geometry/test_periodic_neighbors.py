import itertools

import numpy as np
import pytest

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.geometry import NeighborFinder, PeriodicNeighborGraph
from cristma.structure import CrystalStructure, IndependentSite, SiteComponent
from cristma.symmetry import expand_structure


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def site(label: str, fractional: tuple[float, float, float]) -> IndependentSite:
    return IndependentSite(
        id=f"site:{label}",
        label=label,
        components=(SiteComponent("C", number(1.0)),),
        fractional=tuple(number(value) for value in fractional),
    )


def boundary_crossing_crystal() -> CrystalStructure:
    return CrystalStructure.explicit(
        "boundary",
        UnitCell.cubic(number(10.0)),
        (site("A", (0.05, 0.0, 0.0)), site("B", (0.95, 0.0, 0.0))),
        id="structure:boundary",
    )


def test_periodic_graph_finds_boundary_crossing_image_once() -> None:
    view = expand_structure(boundary_crossing_crystal())

    graph = NeighborFinder(cutoff=1.1).find(view)

    assert isinstance(graph, PeriodicNeighborGraph)
    edge = graph.neighbors(view.atoms[0].id)[0]
    assert edge.target.atom_id == view.atoms[1].id
    assert edge.target.cell_translation == (-1, 0, 0)
    assert edge.distance == pytest.approx(1.0)
    assert edge.vector_cartesian == pytest.approx((-1.0, 0.0, 0.0))
    reverse = graph.neighbors(view.atoms[1].id)[0]
    assert reverse.target.cell_translation == (1, 0, 0)
    assert reverse.vector_cartesian == pytest.approx((1.0, 0.0, 0.0))


def test_periodic_self_neighbors_have_nonzero_translations() -> None:
    crystal = CrystalStructure.explicit(
        "self",
        UnitCell.cubic(number(2.0)),
        (site("A", (0.0, 0.0, 0.0)),),
    )
    view = expand_structure(crystal)

    graph = NeighborFinder(cutoff=2.1).find(view)

    translations = {edge.target.cell_translation for edge in graph.neighbors(view.atoms[0].id)}
    assert (0, 0, 0) not in translations
    assert translations == {
        (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1)
    }


def test_skewed_triclinic_enumeration_is_complete() -> None:
    cell = UnitCell(number(4), number(4), number(4), number(90), number(90), number(20))
    crystal = CrystalStructure.explicit("skew", cell, (site("A", (0, 0, 0)),))
    view = expand_structure(crystal)
    cutoff = 3.0

    graph = NeighborFinder(cutoff=cutoff).find(view)

    expected = {
        translation
        for translation in itertools.product(range(-4, 5), repeat=3)
        if translation != (0, 0, 0)
        and np.linalg.norm(np.asarray(translation) @ cell.matrix) <= cutoff + 1e-12
    }
    actual = {
        edge.target.cell_translation for edge in graph.neighbors(view.atoms[0].id)
    }
    assert (2, -2, 0) in expected
    assert actual == expected


def test_partial_periodicity_never_translates_nonperiodic_axis() -> None:
    crystal = CrystalStructure.explicit(
        "slab",
        UnitCell.cubic(number(10.0)),
        (site("A", (0.05, 0.05, 0.0)), site("B", (0.95, 0.95, 0.0))),
        periodic=(True, True, False),
    )
    view = expand_structure(crystal)

    graph = NeighborFinder(cutoff=1.5).find(view)

    assert graph.neighbors(view.atoms[0].id)
    assert all(
        edge.target.cell_translation[2] == 0
        for atom in view.atoms
        for edge in graph.neighbors(atom.id)
    )
