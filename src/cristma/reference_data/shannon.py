"""Exact Shannon ionic- and crystal-radius reference lookups."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
from functools import lru_cache
from importlib.resources import files
import json

from cristma.chemistry.elements import normalize_element


class ShannonSpinState(StrEnum):
    UNSPECIFIED = "unspecified"
    HIGH_SPIN = "high_spin"
    LOW_SPIN = "low_spin"


@dataclass(frozen=True, slots=True)
class ShannonRadiusRecord:
    symbol: str
    oxidation_state: int
    coordination: str
    spin_state: ShannonSpinState
    crystal_radius: float
    ionic_radius: float
    unit: str = "angstrom"
    dataset_id: str = "cristma.shannon_radii.pymatgen"
    dataset_version: str = "1"


class ShannonRadii:
    """Read-only exact catalog; no CN, oxidation, or spin state is inferred."""

    def __init__(self, records: tuple[ShannonRadiusRecord, ...]) -> None:
        self._records = records

    @classmethod
    @lru_cache(maxsize=1)
    def default(cls) -> "ShannonRadii":
        resource = files("cristma.reference_data").joinpath(
            "resources", "shannon_radii.json"
        )
        payload = json.loads(resource.read_text(encoding="utf-8"))
        records = tuple(
            ShannonRadiusRecord(
                symbol=item["element"],
                oxidation_state=int(item["oxidation_state"]),
                coordination=item["coordination"],
                spin_state=ShannonSpinState(item["spin_state"]),
                crystal_radius=float(item["crystal_radius"]),
                ionic_radius=float(item["ionic_radius"]),
                dataset_id=payload["metadata"]["dataset_id"],
                dataset_version=payload["metadata"]["version"],
            )
            for item in payload["records"]
        )
        return cls(records)

    def find(
        self,
        symbol: str,
        *,
        oxidation_state: int,
        coordination: str | None = None,
        spin_state: ShannonSpinState | str | None = None,
    ) -> tuple[ShannonRadiusRecord, ...]:
        normalized_symbol = normalize_element(symbol)
        normalized_coordination = coordination.strip().upper() if coordination else None
        normalized_spin = ShannonSpinState(spin_state) if spin_state is not None else None
        return tuple(
            record
            for record in self._records
            if record.symbol == normalized_symbol
            and record.oxidation_state == oxidation_state
            and (normalized_coordination is None or record.coordination == normalized_coordination)
            and (normalized_spin is None or record.spin_state is normalized_spin)
        )

    def get_exact(
        self,
        symbol: str,
        *,
        oxidation_state: int,
        coordination: str,
        spin_state: ShannonSpinState | str | None = None,
    ) -> ShannonRadiusRecord:
        matches = self.find(
            symbol,
            oxidation_state=oxidation_state,
            coordination=coordination,
            spin_state=spin_state,
        )
        if len(matches) != 1:
            raise LookupError(
                f"Expected one exact Shannon radius, found {len(matches)} Shannon radii"
            )
        return matches[0]

    @property
    def records(self) -> tuple[ShannonRadiusRecord, ...]:
        return self._records


__all__ = ["ShannonRadii", "ShannonRadiusRecord", "ShannonSpinState"]
