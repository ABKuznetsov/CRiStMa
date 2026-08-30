"""Canonical non-periodic molecular structure records."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import TYPE_CHECKING, Mapping, Protocol, runtime_checkable

from .identity import StructureProvenance
from .occupation import SiteComponent

if TYPE_CHECKING:
    from .view import AtomicView


def _frozen_metadata(values: Mapping[str, object]) -> Mapping[str, object]:
    return MappingProxyType(dict(values))


@runtime_checkable
class Structure(Protocol):
    """Shared minimum contract for periodic and non-periodic structures."""

    name: str
    periodic: tuple[bool, bool, bool]
    provenance: StructureProvenance

    def atomic_view(self, *, expanded: bool = True) -> AtomicView: ...


@dataclass(frozen=True, slots=True)
class MolecularAtom:
    """One finite atom in Cartesian coordinates."""

    id: str
    label: str
    components: tuple[SiteComponent, ...]
    cartesian: tuple[float, float, float]
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        if not self.id or not self.label:
            raise ValueError("atom id and label must not be empty")
        if not self.components:
            raise ValueError("molecular atom must contain at least one component")
        if math.fsum(float(component.occupancy.value) for component in self.components) > 1.0 + 1e-12:
            raise ValueError("molecular atom occupancy exceeds one")
        if len(self.cartesian) != 3 or not all(math.isfinite(value) for value in self.cartesian):
            raise ValueError("atom Cartesian coordinates must contain three finite values")
        object.__setattr__(self, "metadata", _frozen_metadata(self.metadata))


@dataclass(frozen=True, slots=True)
class MolecularBond:
    """A bond connecting two atoms by stable identity."""

    id: str
    atom1_id: str
    atom2_id: str
    order: float | str
    aromatic: bool = False
    stereo: str | None = None

    def __post_init__(self) -> None:
        if not self.id or not self.atom1_id or not self.atom2_id:
            raise ValueError("bond and endpoint IDs must not be empty")


@dataclass(frozen=True, slots=True)
class MolecularGroup:
    """A named atom grouping; grouping alone never implies rigidity."""

    id: str
    label: str
    atom_ids: tuple[str, ...]
    rigid: bool = False
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        if not self.id or not self.label:
            raise ValueError("group id and label must not be empty")
        if not self.atom_ids:
            raise ValueError("molecular group must contain at least one atom")
        if len(set(self.atom_ids)) != len(self.atom_ids):
            raise ValueError("molecular group contains duplicate atom IDs")
        object.__setattr__(self, "metadata", _frozen_metadata(self.metadata))


@dataclass(frozen=True, slots=True)
class MolecularStructure:
    """A finite molecular structure without an artificial periodic cell."""

    name: str
    atoms: tuple[MolecularAtom, ...]
    bonds: tuple[MolecularBond, ...] = ()
    groups: tuple[MolecularGroup, ...] = ()
    id: str | None = None
    periodic: tuple[bool, bool, bool] = (False, False, False)
    provenance: StructureProvenance = field(default_factory=StructureProvenance)
    metadata: Mapping[str, object] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        if any(self.periodic):
            raise ValueError("molecular structure must be non-periodic")
        atom_ids = tuple(atom.id for atom in self.atoms)
        if len(set(atom_ids)) != len(atom_ids):
            raise ValueError("molecular atom IDs must be unique")
        known_atoms = set(atom_ids)

        bond_ids = tuple(bond.id for bond in self.bonds)
        if len(set(bond_ids)) != len(bond_ids):
            raise ValueError("molecular bond IDs must be unique")
        for bond in self.bonds:
            missing = {bond.atom1_id, bond.atom2_id} - known_atoms
            if missing:
                raise ValueError(f"bond references unknown atom: {sorted(missing)!r}")

        group_ids = tuple(group.id for group in self.groups)
        if len(set(group_ids)) != len(group_ids):
            raise ValueError("molecular group IDs must be unique")
        for group in self.groups:
            missing = set(group.atom_ids) - known_atoms
            if missing:
                raise ValueError(f"group references unknown atom: {sorted(missing)!r}")
        object.__setattr__(self, "metadata", _frozen_metadata(self.metadata))

    @property
    def cell(self) -> None:
        return None

    def atomic_view(self, *, expanded: bool = True) -> AtomicView:
        """Expose molecule atoms through the shared numerical structure view."""

        from .properties import AtomicPropertyTable
        from .view import AtomicView

        return AtomicView(
            atoms=self.atoms,
            cell=None,
            periodic=self.periodic,
            properties=AtomicPropertyTable(len(self.atoms)),
        )


__all__ = [
    "MolecularAtom",
    "MolecularBond",
    "MolecularGroup",
    "MolecularStructure",
    "Structure",
]
