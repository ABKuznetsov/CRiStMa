"""Immutable neutral-X-ray structure-factor inputs and results."""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
import math

from cristma.reference_data import (
    NeutralAtomFormFactorTable,
    XRayFormFactorProvenance,
)

from .models import MillerIndex, ReflectionSet, ReflectionSetStatus


def _nonempty(value: str, name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must not be empty")


def _sha256(value: str, name: str) -> None:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")


@dataclass(frozen=True, slots=True)
class XRayScatteringContext:
    """An explicit neutral-atom X-ray scattering model."""

    form_factors: NeutralAtomFormFactorTable = field(repr=False, compare=False)
    model: str = field(default="neutral_atom_f0", init=False)
    anomalous: bool = field(default=False, init=False)

    @classmethod
    @lru_cache(maxsize=1)
    def default(cls) -> "XRayScatteringContext":
        return cls(NeutralAtomFormFactorTable.default())

    @property
    def table_id(self) -> str:
        return self.provenance.dataset_id

    @property
    def table_version(self) -> str:
        return self.provenance.dataset_version

    @property
    def provenance(self) -> XRayFormFactorProvenance:
        return self.form_factors.provenance


@dataclass(frozen=True, slots=True)
class StructureFactorProvenance:
    method: str
    version: str
    space_group_setting_id: int
    cell_fingerprint: str
    table_id: str
    table_version: str
    raw_f_complex: complex
    contribution_scale: float
    extinction_tolerance: float
    normalized_to_zero: bool

    def __post_init__(self) -> None:
        _nonempty(self.method, "structure-factor method")
        _nonempty(self.version, "structure-factor version")
        if (
            isinstance(self.space_group_setting_id, bool)
            or not isinstance(self.space_group_setting_id, int)
            or not 1 <= self.space_group_setting_id <= 530
        ):
            raise ValueError("space-group setting ID must be an integer from 1 to 530")
        _sha256(self.cell_fingerprint, "cell fingerprint")
        _nonempty(self.table_id, "form-factor table ID")
        _nonempty(self.table_version, "form-factor table version")
        if not isinstance(self.raw_f_complex, complex):
            raise TypeError("raw structure factor must be complex")
        if not (
            math.isfinite(self.raw_f_complex.real)
            and math.isfinite(self.raw_f_complex.imag)
        ):
            raise ValueError("raw structure factor must be finite")
        for value, name in (
            (self.contribution_scale, "contribution scale"),
            (self.extinction_tolerance, "extinction tolerance"),
        ):
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError(f"{name} must be a real number")
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{name} must be finite and non-negative")
        if type(self.normalized_to_zero) is not bool:
            raise TypeError("normalized_to_zero must be bool")

    @property
    def raw_amplitude(self) -> float:
        return abs(self.raw_f_complex)


@dataclass(frozen=True, slots=True)
class StructureFactor:
    reflection_id: str
    representative_hkl: MillerIndex
    f_complex: complex
    provenance: StructureFactorProvenance

    def __post_init__(self) -> None:
        _nonempty(self.reflection_id, "reflection ID")
        if not isinstance(self.representative_hkl, MillerIndex):
            raise TypeError("representative_hkl must be MillerIndex")
        if not isinstance(self.f_complex, complex):
            raise TypeError("structure factor must be complex")
        if not (
            math.isfinite(self.f_complex.real) and math.isfinite(self.f_complex.imag)
        ):
            raise ValueError("structure factor must be finite")
        if not isinstance(self.provenance, StructureFactorProvenance):
            raise TypeError("provenance must be StructureFactorProvenance")
        if self.provenance.normalized_to_zero:
            if self.f_complex != 0j:
                raise ValueError("normalized extinct structure factor must be exact zero")
        elif self.f_complex != self.provenance.raw_f_complex:
            raise ValueError("non-normalized structure factor must retain the raw value")

    @property
    def amplitude(self) -> float:
        return abs(self.f_complex)

    @property
    def f_squared(self) -> float:
        return self.amplitude**2


@dataclass(frozen=True, slots=True)
class StructureFactorSet:
    structure_factors: tuple[StructureFactor, ...]
    reflection_set: ReflectionSet
    context: XRayScatteringContext

    def __post_init__(self) -> None:
        if not isinstance(self.reflection_set, ReflectionSet):
            raise TypeError("reflection_set must be ReflectionSet")
        if not isinstance(self.context, XRayScatteringContext):
            raise TypeError("context must be XRayScatteringContext")
        expected = tuple(
            (item.reflection_id, item.representative_hkl)
            for item in self.reflection_set.reflections
        )
        observed = tuple(
            (item.reflection_id, item.representative_hkl)
            for item in self.structure_factors
        )
        if observed != expected:
            raise ValueError(
                "structure factors must exactly match ReflectionSet order and identity"
            )
        for item in self.structure_factors:
            provenance = item.provenance
            if (
                provenance.space_group_setting_id
                != self.reflection_set.space_group_setting_id
                or provenance.cell_fingerprint != self.cell_fingerprint
            ):
                raise ValueError("structure factor belongs to another diffraction input")
            if (
                provenance.table_id != self.context.table_id
                or provenance.table_version != self.context.table_version
            ):
                raise ValueError("structure factor belongs to another scattering context")

    @property
    def space_group_setting_id(self) -> int:
        return self.reflection_set.space_group_setting_id

    @property
    def cell_fingerprint(self) -> str:
        return self.reflection_set.provenance.cell_fingerprint

    @property
    def status(self) -> ReflectionSetStatus:
        return self.reflection_set.status


__all__ = [
    "StructureFactor",
    "StructureFactorProvenance",
    "StructureFactorSet",
    "XRayScatteringContext",
]
