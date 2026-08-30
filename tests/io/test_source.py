import bz2
import gzip
import lzma

import pytest

from cristma.io.source import MappingSourceResolver, decode_bytes, decode_source


def test_gzip_source_is_decoded_with_inner_suffix(tmp_path) -> None:
    path = tmp_path / "POSCAR.gz"
    path.write_bytes(gzip.compress(b"title\n1.0\n"))

    source = decode_source(path)

    assert source.text == "title\n1.0\n"
    assert source.logical_name.endswith("POSCAR")
    assert source.compression == "gzip"


def test_mapping_resolver_is_explicit_and_blocks_parent_escape() -> None:
    resolver = MappingSourceResolver({"POTCAR": b"TITEL = PAW_PBE Si"})

    assert resolver.resolve("POTCAR", from_source="POSCAR").raw.startswith(b"TITEL")
    assert resolver.resolve("../POTCAR", from_source="POSCAR") is None


@pytest.mark.parametrize(
    ("compress", "name", "expected"),
    (
        (gzip.compress, "model.cif.gz", "gzip"),
        (bz2.compress, "model.cif.bz2", "bzip2"),
        (lzma.compress, "model.cif.xz", "xz"),
    ),
)
def test_compression_is_detected_by_magic_bytes(compress, name, expected) -> None:
    source = decode_bytes(compress(b"data_demo\n"), name)

    assert source.text == "data_demo\n"
    assert source.logical_name == "model.cif"
    assert source.compression == expected
