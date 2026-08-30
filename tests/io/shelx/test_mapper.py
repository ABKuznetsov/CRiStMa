import pytest

from cristma.io.shelx.mapper import map_shelx_structures
from cristma.io.shelx.occupancy import ShelxOccupancyExpression
from cristma.io.shelx.parser import parse_shelx


def _map(source: str):
    parsed = parse_shelx(source, source_name="demo.res")
    structures, diagnostics = map_shelx_structures(parsed.document)
    return structures, (*parsed.diagnostics, *diagnostics)


def test_mapper_builds_canonical_structure_with_symbolic_occupancy_context() -> None:
    source = (
        "TITL demo phase in P 1\n"
        "CELL 0.71073 10 11 12 90 91 92\n"
        "LATT -1\n"
        "SFAC Si O\n"
        "FVAR 0.55 0.75\n"
        "RESI 3 LIG\n"
        "PART 2 21\n"
        "Si1 1 0.1 0.2 0.3 11 0.05\n"
        "O1 2 0.2 0.3 0.4 21 0.04\n"
        "HKLF 4\n"
        "END\n"
        "Q1 1 0.5 0.5 0.5 11 0.05 1.2\n"
    )

    structures, diagnostics = _map(source)
    crystal = structures[0]

    assert not [item for item in diagnostics if item.severity.value == "error"]
    assert crystal.name == "demo phase"
    assert crystal.cell.a.value == 10
    assert crystal.cell.beta.value == 91
    assert crystal.space_group.hm_symbol == "P 1"
    assert len(crystal.space_group.operations) == 1
    assert [site.label for site in crystal.sites] == ["Si1", "O1"]
    assert [site.components[0].element for site in crystal.sites] == ["Si", "O"]
    assert [site.total_occupancy for site in crystal.sites] == pytest.approx([1.0, 0.75])
    expression = crystal.sites[1].components[0].metadata["shelx_occupancy"]
    assert isinstance(expression, ShelxOccupancyExpression)
    assert expression.free_variable_index == 2
    assert crystal.sites[0].metadata["shelx_part"] == 2
    assert crystal.sites[0].metadata["shelx_residue"] == {
        "number": 3,
        "class": "LIG",
    }
    assert crystal.provenance.source.format == "shelx"
    assert crystal.provenance.source.source_name == "demo.res"


def test_atoms_after_hklf_and_q_peaks_never_become_sites() -> None:
    source = (
        "CELL 0.71073 10 10 10 90 90 90\nLATT -1\nSFAC C\n"
        "C1 1 0.1 0.2 0.3 11 0.05\nHKLF 4\n"
        "C2 1 0.2 0.3 0.4 11 0.05\nEND\n"
        "Q1 1 0.3 0.4 0.5 11 0.05 1.2\n"
    )

    structures, diagnostics = _map(source)

    assert not [item for item in diagnostics if item.severity.value == "error"]
    assert [site.label for site in structures[0].sites] == ["C1"]


def test_missing_latt_uses_documented_centrosymmetric_p_default() -> None:
    source = (
        "CELL 0.7 10 10 10 90 90 90\nSFAC C\n"
        "C1 1 0.1 0.2 0.3 11 0.05\nEND\n"
    )

    structures, diagnostics = _map(source)

    assert len(structures[0].space_group.operations) == 2
    defaulted = next(item for item in diagnostics if item.code == "shelx.map.latt_defaulted")
    assert defaulted.recovery == "LATT 1"


@pytest.mark.parametrize(
    ("source", "code"),
    [
        ("LATT -1\nSFAC C\nC1 1 0 0 0 11 0.05\nEND\n", "shelx.map.cell_missing"),
        (
            "CELL 0.7 10 10 10 90 90 90\nLATT -1\nSFAC C\nC1 2 0 0 0 11 0.05\nEND\n",
            "shelx.map.sfac_index_invalid",
        ),
        (
            "CELL 0.7 10 10 10 90 90 90\nLATT -1\nSFAC C\nFVAR 0.5\nC1 1 0 0 0 21 0.05\nEND\n",
            "shelx.map.occupancy_invalid",
        ),
    ],
)
def test_invalid_scientific_source_returns_no_partial_structure(
    source: str,
    code: str,
) -> None:
    structures, diagnostics = _map(source)

    assert not structures
    assert code in [item.code for item in diagnostics]
