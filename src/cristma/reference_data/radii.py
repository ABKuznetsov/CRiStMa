"""Exact covalent-radius reference lookups used by structure-search tools."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
import json
from types import MappingProxyType
from typing import Mapping

from cristma.chemistry.elements import normalize_element


@dataclass(frozen=True, slots=True)
class CovalentRadiusRecord:
    symbol: str
    value: float
    variant: str = "unspecified"
    comment: str = ""
    is_default: bool = True
    unit: str = "angstrom"
    dataset_id: str = "cristma.covalent_radii.cordero_2008"
    dataset_version: str = "1"


class CovalentRadii:
    """Exact lookup with no guessed fallback for missing elements."""

    def __init__(self, records: tuple[CovalentRadiusRecord, ...]) -> None:
        self._all_records = records
        default_records = {record.symbol: record for record in records if record.is_default}
        self._records: Mapping[str, CovalentRadiusRecord] = MappingProxyType(
            default_records
        )

    @classmethod
    @lru_cache(maxsize=1)
    def default(cls) -> "CovalentRadii":
        resource = files("cristma.reference_data").joinpath(
            "resources", "covalent_radii.json"
        )
        payload = json.loads(resource.read_text(encoding="utf-8"))
        metadata = payload["metadata"]
        return cls(tuple(
            CovalentRadiusRecord(
                symbol=item["symbol"],
                value=float(item["value"]),
                variant=item["variant"],
                comment=item["comment"],
                is_default=bool(item["default"]),
                dataset_id=metadata["dataset_id"],
                dataset_version=metadata["version"],
            )
            for item in payload["records"]
        ))

    def find(self, symbol: str) -> CovalentRadiusRecord:
        normalized = normalize_element(symbol)
        try:
            return self._records[normalized]
        except KeyError as exc:
            raise KeyError(f"No covalent radius for {normalized}") from exc

    def find_variants(self, symbol: str) -> tuple[CovalentRadiusRecord, ...]:
        normalized = normalize_element(symbol)
        variants = tuple(
            record for record in self._all_records if record.symbol == normalized
        )
        if not variants:
            raise KeyError(f"No covalent radius for {normalized}")
        return variants

    @property
    def records(self) -> tuple[CovalentRadiusRecord, ...]:
        return self._all_records


__all__ = ["CovalentRadii", "CovalentRadiusRecord"]
