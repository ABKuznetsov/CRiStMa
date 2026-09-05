"""Compile the pinned pymatgen Shannon CSV into CrIStMa reference data."""

from __future__ import annotations

import argparse
import csv
from datetime import date
import hashlib
import json
from pathlib import Path


SPIN_STATES = {
    "": "unspecified",
    "High Spin": "high_spin",
    "Low Spin": "low_spin",
}


def compile_csv(source: Path, destination: Path) -> None:
    source_bytes = source.read_bytes()
    with source.open(encoding="utf-8-sig", newline="") as handle:
        rows = tuple(csv.DictReader(handle))
    records = []
    previous_element: str | None = None
    previous_charge: str | None = None
    for row in rows:
        element = row["Element"]
        charge = row["Charge"]
        if not charge and element == previous_element:
            charge = previous_charge or ""
        if not element or not charge:
            raise ValueError("Shannon row lacks element or unambiguous oxidation state")
        records.append({
            "element": element,
            "oxidation_state": int(charge),
            "coordination": row["Coordination"],
            "spin_state": SPIN_STATES[row["Spin State"]],
            "crystal_radius": float(row["Crystal Radius"]),
            "ionic_radius": float(row["Ionic Radius"]),
        })
        previous_element = element
        previous_charge = charge
    records.sort(
        key=lambda item: (
            item["element"],
            item["oxidation_state"],
            item["coordination"],
            item["spin_state"],
        )
    )
    payload = {
        "metadata": {
            "dataset_id": "cristma.shannon_radii.pymatgen",
            "version": "1",
            "original_reference": "R. D. Shannon, Acta Cryst. A32 (1976) 751-767",
            "doi": "10.1107/S0567739476001551",
            "source_repository": "https://github.com/materialsproject/pymatgen",
            "source_commit": "0428f232a569ffe6b16fa030d38ea35a56d70fd6",
            "source_path": "dev_scripts/periodic_table_resources/Shannon_Radii.csv",
            "source_license": "MIT",
            "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
            "generated": date.today().isoformat(),
            "record_count": len(records),
            "provenance_note": "The upstream generator notes that the CSV's earlier data provenance is unknown; values are attributed to the pinned pymatgen artifact and cross-referenced to Shannon (1976).",
        },
        "records": records,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    arguments = parser.parse_args()
    compile_csv(arguments.source, arguments.destination)


if __name__ == "__main__":
    main()
