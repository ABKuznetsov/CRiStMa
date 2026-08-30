from cristma.io.shelx.parser import parse_shelx
from cristma.io.shelx.records import (
    ShelxCellInstruction,
    ShelxEndInstruction,
    ShelxFvarInstruction,
    ShelxHklfInstruction,
    ShelxLattInstruction,
    ShelxPartInstruction,
    ShelxResiInstruction,
    ShelxSfacInstruction,
    ShelxSymmInstruction,
    ShelxZerrInstruction,
)


def test_core_instructions_are_typed_without_losing_source_records() -> None:
    source = (
        "CELL 0.71073 17.6787 4.4989 24.7670 90.000 90.987 90.000\n"
        "ZERR 5.0000 0.0021 0.0005 0.0021 0.000 0.008 0.000\n"
        "LATT 1\n"
        "SYMM 0.5-X, 0.5+Y, 0.5-Z\n"
        "SFAC C H N O\n"
        "FVAR 0.55079 0.25\n"
        "PART 2 21.0\n"
        "RESI 3 LIG\n"
        "HKLF 4\n"
        "END\n"
    )

    result = parse_shelx(source)
    records = result.document.records

    assert result.ok
    assert [type(record) for record in records] == [
        ShelxCellInstruction,
        ShelxZerrInstruction,
        ShelxLattInstruction,
        ShelxSymmInstruction,
        ShelxSfacInstruction,
        ShelxFvarInstruction,
        ShelxPartInstruction,
        ShelxResiInstruction,
        ShelxHklfInstruction,
        ShelxEndInstruction,
    ]
    cell = records[0]
    assert cell.wavelength.value == 0.71073
    assert cell.wavelength.unit == "angstrom"
    assert cell.cell.a.value == 17.6787
    assert cell.cell.beta.value == 90.987
    assert records[1].formula_units.value == 5.0
    assert tuple(value.value for value in records[1].cell_uncertainties) == (
        0.0021,
        0.0005,
        0.0021,
        0.0,
        0.008,
        0.0,
    )
    assert records[2].code == 1
    assert records[4].entries == ("C", "H", "N", "O")
    assert tuple(value.value for value in records[5].values) == (0.55079, 0.25)
    assert (records[6].part, records[6].occupancy_code) == (2, "21.0")
    assert (records[7].residue_number, records[7].residue_class) == (3, "LIG")
    assert records[8].code == 4


def test_invalid_typed_instruction_is_retained_with_diagnostic() -> None:
    result = parse_shelx("CELL invalid\nEND\n")

    assert result.document.records[0].keyword == "CELL"
    assert [item.code for item in result.diagnostics] == ["shelx.parse.invalid_cell"]
    assert not result.ok
