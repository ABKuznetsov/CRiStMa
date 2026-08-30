"""Minimal structural typing contract for finite atomic positions."""

from __future__ import annotations

from typing import Protocol

from .occupation import SiteComponent


class AtomicPosition(Protocol):
    """Coordinate and occupation capabilities shared by atomic row types."""

    id: str
    cartesian: tuple[float, float, float]
    components: tuple[SiteComponent, ...]


__all__ = ["AtomicPosition"]
