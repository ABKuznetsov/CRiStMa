import pytest

from cristma.core.values import MeasuredValue
from cristma.geometry import NeighborFinder, NeighborGraph
from cristma.structure import MolecularAtom, MolecularStructure, SiteComponent


def molecular_view(*rows: tuple[str, tuple[float, float, float]]):
    component = SiteComponent("C", MeasuredValue(1.0, None, "1"))
    atoms = tuple(
        MolecularAtom(atom_id, atom_id, (component,), cartesian)
        for atom_id, cartesian in rows
    )
    return MolecularStructure("finite", atoms).atomic_view()


def test_finite_cutoff_graph_contains_two_directed_edges() -> None:
    view = molecular_view(
        ("A", (0.0, 0.0, 0.0)),
        ("B", (1.0, 0.0, 0.0)),
        ("C", (3.0, 0.0, 0.0)),
    )

    graph = NeighborFinder(cutoff=1.1).find(view)

    assert isinstance(graph, NeighborGraph)
    assert [(edge.target_atom_id, edge.distance) for edge in graph.neighbors("A")] == [
        ("B", pytest.approx(1.0))
    ]
    assert graph.neighbors("B")[0].vector_cartesian == pytest.approx((-1.0, 0.0, 0.0))
    assert graph.neighbors("C") == ()


def test_finite_graph_rejects_nonpositive_cutoff() -> None:
    with pytest.raises(ValueError, match="cutoff must be positive"):
        NeighborFinder(cutoff=0.0)


def test_finder_configuration_is_inspectable_and_clone_is_immutable() -> None:
    finder = NeighborFinder(cutoff=3.0, tolerance=1e-12)

    clone = finder.clone(cutoff=2.5)

    assert finder.get_config() == {"cutoff": 3.0, "tolerance": 1e-12}
    assert clone.get_config() == {"cutoff": 2.5, "tolerance": 1e-12}
    assert finder.cutoff == 3.0


def test_neighbors_are_deterministic_and_zero_length_pairs_are_diagnostics() -> None:
    view = molecular_view(
        ("A", (0.0, 0.0, 0.0)),
        ("C", (1.0, 0.0, 0.0)),
        ("B", (-1.0, 0.0, 0.0)),
        ("D", (0.0, 0.0, 0.0)),
    )

    graph = NeighborFinder(cutoff=1.1).find(view)

    assert [edge.target_atom_id for edge in graph.neighbors("A")] == ["B", "C"]
    assert not graph.neighbors("D") == ()
    assert "geometry.coincident_positions" in {item.code for item in graph.diagnostics}
    assert all(edge.distance > 0 for atom in graph.atoms for edge in graph.neighbors(atom.id))
