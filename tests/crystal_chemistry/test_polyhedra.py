from __future__ import annotations

import math

import pytest

from cristma.chemistry import GrammarOperation, InteractionPriority
from cristma.core.values import MeasuredValue
from cristma.crystallography import GeometricContact
from cristma.crystal_chemistry import (
    ComponentPairInterpretation,
    ContactClassification,
    CoordinationShell,
    EvidenceStatus,
    PolyhedronBuilder,
    ResolutionStatus,
    ResolvedContact,
)
from cristma.structure import MolecularAtom, MolecularStructure, SiteComponent


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def build_shell(
    vertices: tuple[tuple[float, float, float], ...],
    *,
    status: ResolutionStatus = ResolutionStatus.RESOLVED,
):
    center_component = SiteComponent("Ca", number(1.0))
    ligand_component = SiteComponent("O", number(1.0))
    atoms = [MolecularAtom("center", "Ca", (center_component,), (0.0, 0.0, 0.0))]
    contacts = []
    for index, vertex in enumerate(vertices):
        atom_id = f"ligand:{index}"
        atoms.append(MolecularAtom(atom_id, f"O{index}", (ligand_component,), vertex))
        distance = math.sqrt(sum(value * value for value in vertex))
        geometric = GeometricContact(
            f"contact:{index}", "center", atom_id, None, distance, vertex,
            "site:center", f"site:ligand:{index}", "test",
        )
        interpretation = ComponentPairInterpretation(
            center_component.species, ligand_component.species,
            1.0, 1.0, 2.42, distance / 2.42, 1.0,
            GrammarOperation.CENTRE_LIGAND_SHELL,
            InteractionPriority.PRIMARY,
            ("Ca",), ("O",),
        )
        contacts.append(ResolvedContact(
            geometric,
            GrammarOperation.CENTRE_LIGAND_SHELL,
            InteractionPriority.PRIMARY,
            ContactClassification.PRIMARY,
            (interpretation,),
            interpretation.normalized_distance,
            interpretation.normalized_distance,
            1.0,
            (),
            (),
        ))
    view = MolecularStructure("fixture", tuple(atoms)).atomic_view()
    shell = CoordinationShell(
        "site:center", "center",
        tuple(contacts) if status is ResolutionStatus.RESOLVED else (),
        len(contacts) if status is ResolutionStatus.RESOLVED else 0,
        float(len(contacts)) if status is ResolutionStatus.RESOLVED else 0.0,
        status,
    )
    return shell, view


def test_tetrahedral_shell_builds_closed_hull() -> None:
    shell, view = build_shell(((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)))

    result = PolyhedronBuilder().build(shell, view)

    assert result.status is ResolutionStatus.RESOLVED
    assert result.polyhedron is not None
    assert len(result.polyhedron.faces) == 4
    assert result.polyhedron.volume == pytest.approx(8.0 / 3.0)
    assert result.polyhedron.center_offset == pytest.approx(0.0)


def test_polyhedron_builder_configuration_is_explicit_and_cloneable() -> None:
    builder = PolyhedronBuilder(tolerance=1e-8)

    assert builder.get_config() == {"tolerance": 1e-8}
    assert builder.clone(tolerance=1e-7).tolerance == pytest.approx(1e-7)
    assert builder.tolerance == pytest.approx(1e-8)


def test_octahedral_shell_has_eight_triangular_faces() -> None:
    shell, view = build_shell(((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)))

    polyhedron = PolyhedronBuilder().build(shell, view).polyhedron

    assert polyhedron is not None
    assert len(polyhedron.faces) == 8
    assert {len(face) for face in polyhedron.faces} == {3}


def test_cube_exposes_six_quadrilateral_faces() -> None:
    vertices = tuple(
        (float(x), float(y), float(z))
        for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)
    )
    shell, view = build_shell(vertices)

    polyhedron = PolyhedronBuilder().build(shell, view).polyhedron

    assert polyhedron is not None
    assert len(polyhedron.faces) == 6
    assert {len(face) for face in polyhedron.faces} == {4}
    assert polyhedron.volume == pytest.approx(8.0)


@pytest.mark.parametrize(
    ("vertices", "status", "expected"),
    [
        (((-1, 0, 0), (1, 0, 0)), ResolutionStatus.RESOLVED, ResolutionStatus.NOT_APPLICABLE),
        (((-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0)), ResolutionStatus.RESOLVED, ResolutionStatus.NOT_APPLICABLE),
        (((1, 0, 0),), ResolutionStatus.AMBIGUOUS, ResolutionStatus.AMBIGUOUS),
        (((1, 0, 0),), ResolutionStatus.INCOMPLETE, ResolutionStatus.INCOMPLETE),
    ],
)
def test_non_success_reason_is_preserved(vertices, status, expected) -> None:
    shell, view = build_shell(vertices, status=status)

    result = PolyhedronBuilder().build(shell, view)

    assert result.polyhedron is None
    assert result.status is expected
    assert result.diagnostics


def test_inconsistent_equivalent_topology_returns_diagnostic() -> None:
    tetra_shell, tetra_view = build_shell(((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)))
    octa_shell, octa_view = build_shell(((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)))
    tetrahedron = PolyhedronBuilder().build(tetra_shell, tetra_view).polyhedron
    octahedron = PolyhedronBuilder().build(octa_shell, octa_view).polyhedron
    assert tetrahedron is not None and octahedron is not None

    diagnostics = PolyhedronBuilder().validate_orbit((tetrahedron, octahedron))

    assert {item.code for item in diagnostics} == {
        "crystal_chemistry.polyhedron.symmetry_inconsistent"
    }
