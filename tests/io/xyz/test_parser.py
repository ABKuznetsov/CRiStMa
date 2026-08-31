import numpy as np
import pytest

from cristma.io.xyz import index_xyz
from cristma.io.xyz.parser import load_xyz_frame, validate_xyz_frame
from cristma.structure import FrameReference, SourceReference


def _load(source: str, name: str = "sample.xyz"):
    document, diagnostics = index_xyz(source, name)
    assert not diagnostics
    span = document.frames[0]
    reference = FrameReference(
        index=0,
        source=SourceReference(name, "xyz", "frame:0", span.start_offset, span.end_offset),
    )
    return document, span, reference


def test_all_declared_property_types_are_typed() -> None:
    source = '''2
Properties=species:S:1:Z:I:1:pos:R:3:forces:R:3:fixed:L:1 label="typed"
Si 14 0 0 0 1.0 2.0 3.0 T
O 8 1 0 0 -1.0 -2.0 -3.0 F
'''
    document, _, reference = _load(source)

    frame = load_xyz_frame(document, reference)

    assert frame.columns["species"].dtype.kind in {"U", "O"}
    assert frame.columns["Z"].dtype.kind == "i"
    assert frame.columns["forces"].shape == (2, 3)
    assert frame.columns["fixed"].dtype.kind == "b"
    assert frame.metadata["label"] == "typed"
    assert frame.source is reference.source


def test_plain_xyz_uses_implicit_columns_and_ignores_trailing_tokens() -> None:
    document, span, reference = _load("1\nwater\nO 0 1 2 charge=ignored\n")

    frame = load_xyz_frame(document, reference)

    assert tuple(item.name for item in frame.schema) == ("species", "pos")
    assert frame.columns["species"].tolist() == ["O"]
    assert frame.columns["pos"].tolist() == [[0.0, 1.0, 2.0]]
    assert any(
        item.code == "xyz.map.uninterpreted_plain_columns"
        for item in validate_xyz_frame(document, span)
    )


def test_row_width_must_equal_schema_width() -> None:
    document, _, reference = _load(
        "1\nProperties=species:S:1:pos:R:3\nSi 0 0\n"
    )

    with pytest.raises(ValueError, match="column count"):
        load_xyz_frame(document, reference)


@pytest.mark.parametrize(
    "kind, value",
    [("I", "1.2"), ("R", "not-a-number"), ("L", "maybe")],
)
def test_invalid_typed_cells_are_rejected(kind: str, value: str) -> None:
    document, _, reference = _load(
        f"1\nProperties=species:S:1:value:{kind}:1\nSi {value}\n"
    )

    with pytest.raises(ValueError, match="value"):
        load_xyz_frame(document, reference)


def test_lattice_and_pbc_are_kept_as_explicit_independent_metadata() -> None:
    document, span, reference = _load(
        '1\nLattice="2 0 0 0 2 0 0 0 2" Properties=species:S:1:pos:R:3\nSi 0 0 0\n'
    )

    frame = load_xyz_frame(document, reference)

    assert np.array_equal(frame.lattice, np.diag([2.0, 2.0, 2.0]))
    assert frame.pbc is None
    assert any(item.code == "xyz.map.lattice_without_pbc" for item in validate_xyz_frame(document, span))


def test_validation_reports_schema_or_cell_errors_without_materializing_frame() -> None:
    document, span, _ = _load(
        "1\nProperties=species:S:1:pos:R:3\nSi bad 0 0\n"
    )

    diagnostics = validate_xyz_frame(document, span)

    assert any(item.severity.value == "error" and item.code == "xyz.map.value_invalid" for item in diagnostics)


def test_validation_reports_unknown_species_and_species_z_conflict() -> None:
    unknown_document, unknown_span, _ = _load("1\nunknown\nXx 0 0 0\n")
    conflict_document, conflict_span, _ = _load(
        "1\nProperties=species:S:1:Z:I:1:pos:R:3\nSi 8 0 0 0\n"
    )

    assert any(
        item.code == "xyz.map.species_unresolved"
        for item in validate_xyz_frame(unknown_document, unknown_span)
    )
    assert any(
        item.severity.value == "error" and item.code == "xyz.map.species_conflict"
        for item in validate_xyz_frame(conflict_document, conflict_span)
    )


@pytest.mark.parametrize(
    "schema",
    [
        "forces:R:3",
        "species:S:1:forces:R:3",
        "species:S:1:pos:I:3",
    ],
)
def test_extended_schema_requires_identity_and_cartesian_position(schema: str) -> None:
    width = sum(int(value) for value in schema.split(":")[2::3])
    document, span, _ = _load(
        f"1\nProperties={schema}\n" + " ".join(["0"] * width) + "\n"
    )

    assert any(
        item.severity.value == "error" and item.code == "xyz.map.schema_invalid"
        for item in validate_xyz_frame(document, span)
    )


def test_true_pbc_without_lattice_is_invalid() -> None:
    document, span, _ = _load(
        '1\nProperties=species:S:1:pos:R:3 pbc="T F F"\nSi 0 0 0\n'
    )

    assert any(
        item.severity.value == "error" and item.code == "xyz.map.lattice_required"
        for item in validate_xyz_frame(document, span)
    )
