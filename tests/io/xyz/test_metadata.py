import numpy as np
import pytest

from cristma.io.xyz.metadata import parse_property_schema, parse_xyz_metadata


def test_extxyz_metadata_parses_special_and_unknown_values() -> None:
    metadata = parse_xyz_metadata(
        'Lattice="2 0 0 0 2 0 0 0 2" '
        'Properties=species:S:1:pos:R:3:forces:R:3 '
        'pbc="T T F" energy=-1.25 count=4 label="relaxed cell"'
    )

    assert np.array_equal(metadata.lattice, np.diag([2.0, 2.0, 2.0]))
    assert metadata.pbc == (True, True, False)
    assert metadata.values["energy"] == -1.25
    assert metadata.values["count"] == 4
    assert metadata.values["label"] == "relaxed cell"
    assert [item.name for item in metadata.schema] == ["species", "pos", "forces"]


def test_metadata_supports_escaped_quotes_and_logicals() -> None:
    metadata = parse_xyz_metadata('label="a \\"quoted\\" cell" active=T flags="T F T"')

    assert metadata.values["label"] == 'a "quoted" cell'
    assert metadata.values["active"] is True
    assert metadata.values["flags"] == (True, False, True)


def test_schema_rejects_duplicate_names_or_wrong_triplets() -> None:
    with pytest.raises(ValueError, match="duplicate"):
        parse_property_schema("species:S:1:pos:R:3:pos:R:3")
    with pytest.raises(ValueError, match="triplets"):
        parse_property_schema("species:S:1:broken")


@pytest.mark.parametrize(
    "comment, message",
    [
        ('Lattice="1 0 0"', "Lattice"),
        ('pbc="T F"', "pbc"),
        ('Properties=species:S:1 Properties=pos:R:3', "duplicate"),
    ],
)
def test_malformed_special_metadata_is_rejected(comment: str, message: str) -> None:
    with pytest.raises(ValueError, match=message):
        parse_xyz_metadata(comment)


def test_arbitrary_plain_comment_is_not_mistaken_for_metadata() -> None:
    metadata = parse_xyz_metadata("water optimized at room temperature")

    assert metadata.values == {"comment": "water optimized at room temperature"}
    assert metadata.schema == ()

