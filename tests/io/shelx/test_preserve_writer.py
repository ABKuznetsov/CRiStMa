import pytest

from cristma.io.shelx.parser import parse_shelx
from cristma.io.shelx.writer import write_shelx_document


@pytest.mark.parametrize(
    "source",
    [
        "TITL demo\nREM keep  two spaces\n\nEND\nQ1 1 0 0 0 11 0.05\n",
        "TITL demo\r\nREM keep  two spaces\r\n\r\nEND\r\nQ1 1 0 0 0 11 0.05\r\n",
    ],
)
def test_unchanged_document_round_trips_exactly(source: str) -> None:
    document = parse_shelx(source, source_name="demo.res").document

    assert write_shelx_document(document, mode="preserve") == source


def test_document_writer_rejects_non_preserve_mode() -> None:
    document = parse_shelx("END\n").document

    with pytest.raises(ValueError, match="preserve"):
        write_shelx_document(document, mode="canonical")
