from cristma.io.formats import builtin_format_descriptors
from cristma.io.registry import FormatRegistry


def test_builtin_cif_descriptor_is_content_aware_and_multiple() -> None:
    registry = FormatRegistry(builtin_format_descriptors())

    descriptor = registry.select("data_demo\n_tag value\n")

    assert descriptor.name == "cif"
    assert descriptor.capabilities.multiple
