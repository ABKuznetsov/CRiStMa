"""Declarative, lazily instantiated structure format descriptions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Protocol, runtime_checkable


@runtime_checkable
class FormatHandler(Protocol):
    """Reader implementation created only after its descriptor is selected."""

    name: str
    suffixes: tuple[str, ...]

    def probe(self, source: str) -> float: ...

    def read_text(self, source: str, source_name: str | None = None) -> object: ...


@dataclass(frozen=True, slots=True)
class FormatCapabilities:
    text: bool = True
    binary: bool = False
    multiple: bool = False
    lazy_frames: bool = False


@dataclass(frozen=True, slots=True)
class FormatDescriptor:
    name: str
    aliases: tuple[str, ...]
    suffixes: tuple[str, ...]
    basenames: tuple[str, ...]
    probe: Callable[[str], float]
    factory: Callable[[], FormatHandler]
    capabilities: FormatCapabilities

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("format name must not be empty")


def descriptor_for(handler: FormatHandler) -> FormatDescriptor:
    """Wrap a legacy eager handler in the descriptor contract."""

    return FormatDescriptor(
        name=handler.name,
        aliases=(),
        suffixes=tuple(handler.suffixes),
        basenames=(),
        probe=handler.probe,
        factory=lambda: handler,
        capabilities=FormatCapabilities(text=True),
    )


def _cif_handler() -> FormatHandler:
    from .cif.handler import CifFormatHandler

    return CifFormatHandler()


def _shelx_handler() -> FormatHandler:
    from .shelx.handler import ShelxFormatHandler

    return ShelxFormatHandler()


def _vasp_handler() -> FormatHandler:
    from .vasp.handler import VaspFormatHandler

    return VaspFormatHandler()


def _xyz_handler() -> FormatHandler:
    from .xyz.handler import XyzFormatHandler

    return XyzFormatHandler()


def builtin_format_descriptors() -> tuple[FormatDescriptor, ...]:
    """Return built-ins without importing their parser or mapper modules."""

    from .cif.probe import probe_cif
    from .shelx.probe import probe_shelx
    from .vasp.probe import probe_vasp
    from .xyz.probe import probe_xyz

    return (
        FormatDescriptor(
            name="cif",
            aliases=("mmcif",),
            suffixes=(".cif", ".mmcif", ".mcif"),
            basenames=(),
            probe=probe_cif,
            factory=_cif_handler,
            capabilities=FormatCapabilities(text=True, multiple=True),
        ),
        FormatDescriptor(
            name="shelx",
            aliases=("res", "ins"),
            suffixes=(".res", ".ins"),
            basenames=(),
            probe=probe_shelx,
            factory=_shelx_handler,
            capabilities=FormatCapabilities(text=True, multiple=False),
        ),
        FormatDescriptor(
            name="vasp",
            aliases=("poscar", "contcar", "xdatcar", "outcar", "vasprun"),
            suffixes=(".xml",),
            basenames=("POSCAR", "CONTCAR", "XDATCAR", "OUTCAR", "vasprun.xml"),
            probe=probe_vasp,
            factory=_vasp_handler,
            capabilities=FormatCapabilities(text=True, multiple=True, lazy_frames=True),
        ),
        FormatDescriptor(
            name="xyz",
            aliases=("extxyz",),
            suffixes=(".xyz", ".extxyz"),
            basenames=(),
            probe=probe_xyz,
            factory=_xyz_handler,
            capabilities=FormatCapabilities(text=True, multiple=True, lazy_frames=True),
        ),
    )


__all__ = [
    "FormatCapabilities",
    "FormatDescriptor",
    "FormatHandler",
    "builtin_format_descriptors",
    "descriptor_for",
]
