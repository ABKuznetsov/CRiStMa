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


def builtin_format_descriptors() -> tuple[FormatDescriptor, ...]:
    """Return built-ins without importing their parser or mapper modules."""

    from .cif.probe import probe_cif
    from .shelx.probe import probe_shelx

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
    )


__all__ = [
    "FormatCapabilities",
    "FormatDescriptor",
    "FormatHandler",
    "builtin_format_descriptors",
    "descriptor_for",
]
