import numpy as np
import pytest

from cristma.structure import SourceReference
from cristma.io.xyz.document import (
    XyzDocument,
    XyzFrame,
    XyzFrameSpan,
    XyzPropertySpec,
)


def test_property_spec_requires_known_type_and_positive_width() -> None:
    assert XyzPropertySpec("forces", "R", 3).width == 3
    with pytest.raises(ValueError, match="name"):
        XyzPropertySpec("", "R", 3)
    with pytest.raises(ValueError, match="type"):
        XyzPropertySpec("forces", "Q", 3)
    with pytest.raises(ValueError, match="width"):
        XyzPropertySpec("forces", "R", 0)


def test_frame_arrays_and_metadata_are_immutable() -> None:
    frame = XyzFrame(
        name="water",
        atom_count=1,
        comment="demo",
        metadata={"energy": -1.0},
        schema=(XyzPropertySpec("pos", "R", 3),),
        columns={"pos": np.array([[0.0, 0.0, 0.0]])},
        lattice=None,
        pbc=None,
        source=SourceReference("water.xyz", "xyz", "frame:0", 0, 20),
    )

    with pytest.raises(ValueError):
        frame.columns["pos"][0, 0] = 1.0
    with pytest.raises(TypeError):
        frame.metadata["energy"] = 0.0


def test_frame_rejects_property_row_count_mismatch() -> None:
    with pytest.raises(ValueError, match="atom count"):
        XyzFrame(
            name="bad",
            atom_count=2,
            comment="",
            metadata={},
            schema=(XyzPropertySpec("pos", "R", 3),),
            columns={"pos": np.array([[0.0, 0.0, 0.0]])},
            lattice=None,
            pbc=None,
            source=SourceReference("bad.xyz", "xyz", "frame:0", 0, 10),
        )


def test_document_validates_span_order() -> None:
    with pytest.raises(ValueError, match="ordered"):
        XyzDocument(
            "source",
            "bad.xyz",
            (
                XyzFrameSpan(1, 0, 2, 3, 4, 5, 6),
                XyzFrameSpan(0, 0, 0, 1, 1, 2, 2),
            ),
        )
