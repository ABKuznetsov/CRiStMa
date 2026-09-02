from __future__ import annotations


def test_crystal_chemistry_tools_import_from_package() -> None:
    from cristma.crystallography import GeometricContact, geometric_contacts
    from cristma.crystal_chemistry import (
        CoordinationShellResolver,
        PolyhedronBuilder,
        ShellResolutionPolicy,
    )

    assert GeometricContact is not None
    assert geometric_contacts is not None
    assert CoordinationShellResolver is not None
    assert PolyhedronBuilder is not None
    assert ShellResolutionPolicy is not None
