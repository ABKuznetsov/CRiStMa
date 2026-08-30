import pytest

from cristma.io.shelx.mapper import map_shelx_structures
from cristma.io.shelx.parser import parse_shelx


def test_u_iso_and_shelx_u_aniso_order_map_to_canonical_tensor() -> None:
    source = (
        "CELL 0.71073 10 10 10 90 90 90\nLATT -1\nSFAC C O\n"
        "C1 1 0.1 0.2 0.3 11 0.05\n"
        "O1 2 0.2 0.3 0.4 11 0.01 0.02 0.03 0.004 0.005 0.006\n"
        "END\n"
    )

    structures, diagnostics = map_shelx_structures(parse_shelx(source).document)
    c_site, o_site = structures[0].sites

    assert not [item for item in diagnostics if item.severity.value == "error"]
    assert c_site.displacement.kind == "U_iso"
    assert c_site.displacement.isotropic.value == 0.05
    tensor = o_site.displacement.tensor
    assert [[value.value for value in row] for row in tensor] == [
        [0.01, 0.006, 0.005],
        [0.006, 0.02, 0.004],
        [0.005, 0.004, 0.03],
    ]


def test_anisotropic_tensor_is_rotated_for_symmetry_image() -> None:
    source = (
        "CELL 0.71073 10 10 10 90 90 90\nLATT -1\nSYMM y,x,z\nSFAC O\n"
        "O1 1 0.1 0.2 0.3 11 0.01 0.02 0.03 0 0 0\nEND\n"
    )

    structures, diagnostics = map_shelx_structures(parse_shelx(source).document)
    atoms = structures[0].atomic_view().atoms

    assert not [item for item in diagnostics if item.severity.value == "error"]
    assert len(atoms) == 2
    assert atoms[0].displacement.tensor[0][0].value == pytest.approx(0.01)
    assert atoms[1].displacement.tensor[0][0].value == pytest.approx(0.02)
    assert atoms[1].displacement.tensor[1][1].value == pytest.approx(0.01)


def test_negative_u_dependency_is_not_misreported_as_physical_displacement() -> None:
    source = (
        "CELL 0.71073 10 10 10 90 90 90\nLATT -1\nSFAC H\n"
        "H1 1 0.1 0.2 0.3 11 -1.2\nEND\n"
    )

    structures, diagnostics = map_shelx_structures(parse_shelx(source).document)

    assert structures[0].sites[0].displacement is None
    assert "shelx.map.displacement_dependency_unmapped" in [
        item.code for item in diagnostics
    ]
