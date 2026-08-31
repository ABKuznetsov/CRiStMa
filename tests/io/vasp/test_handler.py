from cristma.io.formats import builtin_format_descriptors
from cristma.io.registry import FormatRegistry


POSCAR_TEXT = """Silicon
1.0
5.43 0 0
0 5.43 0
0 0 5.43
Si
1
Direct
0 0 0
"""


def test_poscar_content_selects_lazy_vasp_descriptor_without_filename() -> None:
    registry = FormatRegistry(builtin_format_descriptors())

    descriptor = registry.select(POSCAR_TEXT, basename="renamed.data")

    assert descriptor.name == "vasp"
    assert descriptor.capabilities.multiple
    assert descriptor.capabilities.lazy_frames


def test_handler_maps_poscar_to_canonical_structure() -> None:
    registry = FormatRegistry(builtin_format_descriptors())

    result = registry.read_text(POSCAR_TEXT, source_name="POSCAR")

    assert result.ok
    assert result.structures[0].name == "Silicon"
    assert result.structures[0].sites[0].components[0].element == "Si"
