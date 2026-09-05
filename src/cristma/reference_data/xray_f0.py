"""Immutable neutral-atom X-ray form factors on the CrIStMa s convention."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import hashlib
from importlib.resources import files
import json
import math
from numbers import Integral, Real
from types import MappingProxyType
from typing import Mapping

import numpy as np


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")


@dataclass(frozen=True, slots=True)
class XRayFormFactorProvenance:
    dataset_id: str
    dataset_version: str
    source: str
    xraylib_version: str
    xraylib_commit: str
    source_sha256: str
    data_sha256: str
    q_convention: str
    cristma_variable: str
    s_unit: str
    interpolation: str
    supported_atomic_numbers: tuple[int, int]
    generator: str


@dataclass(frozen=True, slots=True)
class _FormFactorSeries:
    s: np.ndarray
    f0: np.ndarray
    second_derivatives: np.ndarray


class NeutralAtomFormFactorTable:
    """Read-only xraylib FF_Rayl table evaluated as f0(s)."""

    def __init__(
        self,
        provenance: XRayFormFactorProvenance,
        series: Mapping[int, _FormFactorSeries],
    ) -> None:
        self._provenance = provenance
        self._series = MappingProxyType(dict(series))

    @classmethod
    @lru_cache(maxsize=1)
    def default(cls) -> "NeutralAtomFormFactorTable":
        resource = files("cristma.reference_data").joinpath(
            "resources", "xray", "xray_f0.json"
        )
        payload = json.loads(resource.read_text(encoding="ascii"))
        metadata = payload["metadata"]
        elements = payload["elements"]
        data_sha256 = hashlib.sha256(_canonical_bytes(elements)).hexdigest()
        if data_sha256 != metadata["data_sha256"]:
            raise ValueError("X-ray form-factor resource checksum mismatch")
        supported = tuple(metadata["supported_atomic_numbers"])
        if supported != (1, 98):
            raise ValueError("X-ray form-factor resource has an unsupported Z range")

        series: dict[int, _FormFactorSeries] = {}
        for item in elements:
            atomic_number = int(item["atomic_number"])
            arrays = tuple(
                np.asarray(item[name], dtype=float)
                for name in ("s", "f0", "second_derivatives")
            )
            if not arrays[0].size or len({array.size for array in arrays}) != 1:
                raise ValueError("X-ray form-factor arrays have inconsistent lengths")
            if not all(np.all(np.isfinite(array)) for array in arrays):
                raise ValueError("X-ray form-factor resource contains non-finite values")
            if not np.all(np.diff(arrays[0]) > 0):
                raise ValueError("X-ray form-factor s grid must be strictly increasing")
            for array in arrays:
                array.setflags(write=False)
            series[atomic_number] = _FormFactorSeries(*arrays)
        if tuple(sorted(series)) != tuple(range(1, 99)):
            raise ValueError("X-ray form-factor resource must contain Z=1..98")

        provenance = XRayFormFactorProvenance(
            dataset_id=metadata["dataset_id"],
            dataset_version=metadata["version"],
            source=metadata["source"],
            xraylib_version=metadata["xraylib_version"],
            xraylib_commit=metadata["xraylib_commit"],
            source_sha256=metadata["source_sha256"],
            data_sha256=metadata["data_sha256"],
            q_convention=metadata["q_convention"],
            cristma_variable=metadata["cristma_variable"],
            s_unit=metadata["s_unit"],
            interpolation=metadata["interpolation"],
            supported_atomic_numbers=supported,
            generator=metadata["generator"],
        )
        return cls(provenance, series)

    load_default = default

    @property
    def provenance(self) -> XRayFormFactorProvenance:
        return self._provenance

    def evaluate(self, atomic_number: int, s: float) -> float:
        if isinstance(atomic_number, bool) or not isinstance(atomic_number, Integral):
            raise TypeError("atomic number must be an integer")
        normalized_z = int(atomic_number)
        if normalized_z not in self._series:
            raise ValueError("atomic number is outside the supported range 1..98")
        if isinstance(s, bool) or not isinstance(s, Real):
            raise TypeError("s must be a real number")
        normalized_s = float(s)
        if not math.isfinite(normalized_s) or normalized_s < 0:
            raise ValueError("s must be finite and non-negative")
        if normalized_s == 0.0:
            return float(normalized_z)

        series = self._series[normalized_z]
        if normalized_s < series.s[0] or normalized_s > series.s[-1] + 1e-7:
            raise ValueError("s is outside the tabulated range")
        normalized_s = min(normalized_s, float(series.s[-1]))
        upper = int(np.searchsorted(series.s, normalized_s, side="right"))
        upper = min(max(upper, 1), series.s.size - 1)
        lower = upper - 1
        width = float(series.s[upper] - series.s[lower])
        if width == 0.0:
            return float((series.f0[lower] + series.f0[upper]) / 2.0)
        a = float((series.s[upper] - normalized_s) / width)
        b = float((normalized_s - series.s[lower]) / width)
        value = (
            a * series.f0[lower]
            + b * series.f0[upper]
            + (
                (a**3 - a) * series.second_derivatives[lower]
                + (b**3 - b) * series.second_derivatives[upper]
            )
            * width**2
            / 6.0
        )
        return float(value)


__all__ = ["NeutralAtomFormFactorTable", "XRayFormFactorProvenance"]
