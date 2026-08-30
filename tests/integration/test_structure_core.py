from pathlib import Path

import cristma
from cristma.geometry import CoordinationAnalyzer, NeighborFinder
from cristma.structure import AtomicPosition, AtomicView, ExpandedAtom, PeriodicAtomRef
from cristma.symmetry import SymmetryImageProvenance, expand_structure


FIXTURES = Path("tests/fixtures/cif")


def test_real_cif_traverses_structure_core() -> None:
    result = cristma.read(FIXTURES / "cod_3000098_barium_borate.cif")

    assert result.ok
    crystal = result.structures[0]
    view = expand_structure(crystal)
    graph = NeighborFinder(cutoff=2.6).find(view)
    coordination = CoordinationAnalyzer().analyze(view, graph)

    assert view.atoms
    assert len(graph.atoms) == len(view.atoms)
    assert len(coordination.environments) == len(view.atoms)
    assert all(atom.source_site_id for atom in view.atoms)


def test_mixed_disorder_is_one_position_without_zero_length_edges() -> None:
    crystal = cristma.read(FIXTURES / "mixed_disorder.cif").structures[0]

    view = expand_structure(crystal)
    graph = NeighborFinder(cutoff=2.6).find(view)

    assert len(view.atoms) == 1
    assert tuple(component.element for component in view.atoms[0].components) == ("La", "Zr")
    assert all(edge.distance > 0 for edge in graph.neighbors(view.atoms[0].id))


def test_structure_core_namespace_contract_is_importable() -> None:
    assert AtomicPosition and AtomicView and ExpandedAtom and PeriodicAtomRef
    assert SymmetryImageProvenance and expand_structure
