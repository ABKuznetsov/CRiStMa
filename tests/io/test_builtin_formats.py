import os
from pathlib import Path
import subprocess
import sys

from cristma.io.formats import builtin_format_descriptors
from cristma.io.registry import FormatRegistry


def test_builtin_cif_descriptor_is_content_aware_and_multiple() -> None:
    registry = FormatRegistry(builtin_format_descriptors())

    descriptor = registry.select("data_demo\n_tag value\n")

    assert descriptor.name == "cif"
    assert descriptor.capabilities.multiple


def test_builtin_shelx_descriptor_is_content_aware_and_single_structure() -> None:
    registry = FormatRegistry(builtin_format_descriptors())

    descriptor = registry.select(
        "TITL demo\nCELL 0.71073 10 10 10 90 90 90\nLATT -1\nSFAC C\nEND\n"
    )

    assert descriptor.name == "shelx"
    assert descriptor.aliases == ("res", "ins")
    assert descriptor.suffixes == (".res", ".ins")
    assert not descriptor.capabilities.multiple


def test_shelx_suffix_selects_descriptor_for_partial_source() -> None:
    registry = FormatRegistry(builtin_format_descriptors())

    assert registry.select("TITL partial\n", suffix=".res").name == "shelx"
    assert registry.select("TITL partial\n", suffix=".ins").name == "shelx"


def test_builtin_descriptor_does_not_import_shelx_parser_or_mapper() -> None:
    root = Path(__file__).parents[2]
    script = (
        "import sys; "
        "from cristma.io.formats import builtin_format_descriptors; "
        "builtin_format_descriptors(); "
        "assert 'cristma.io.shelx.parser' not in sys.modules; "
        "assert 'cristma.io.shelx.mapper' not in sys.modules"
    )
    environment = dict(os.environ, PYTHONPATH=str(root / "src"))

    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr


def test_builtin_descriptor_does_not_import_vasp_implementations() -> None:
    root = Path(__file__).parents[2]
    script = (
        "import sys; "
        "from cristma.io.formats import builtin_format_descriptors; "
        "builtin_format_descriptors(); "
        "assert 'cristma.io.vasp.poscar' not in sys.modules; "
        "assert 'cristma.io.vasp.outcar' not in sys.modules; "
        "assert 'cristma.io.vasp.vasprun' not in sys.modules"
    )
    environment = dict(os.environ, PYTHONPATH=str(root / "src"))

    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr


def test_xyz_descriptor_is_lazy_and_content_aware() -> None:
    registry = FormatRegistry(builtin_format_descriptors())

    descriptor = registry.select("1\nwater\nO 0 0 0\n")

    assert descriptor.name == "xyz"
    assert descriptor.aliases == ("extxyz",)
    assert descriptor.capabilities.multiple
    assert descriptor.capabilities.lazy_frames


def test_builtin_descriptor_does_not_import_xyz_parser_or_mapper() -> None:
    root = Path(__file__).parents[2]
    script = (
        "import sys; "
        "from cristma.io.formats import builtin_format_descriptors; "
        "builtin_format_descriptors(); "
        "assert 'cristma.io.xyz.parser' not in sys.modules; "
        "assert 'cristma.io.xyz.mapper' not in sys.modules"
    )
    environment = dict(os.environ, PYTHONPATH=str(root / "src"))

    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=root,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
