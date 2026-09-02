"""Compile pinned QCElemental/Cordero covalent radii for CRiStMa."""

from __future__ import annotations

import argparse
import ast
from datetime import date
import hashlib
import json
from pathlib import Path


SOURCE_SHA256 = "9ac22bedfc04ead3567ebf0484fe09583e959d679aec99b69a0aef13388cb63e"
SOURCE_COMMIT = "c4eb31cff9c7041f4767804a0076e35343df8177"


def _source_payload(source: Path) -> dict[str, object]:
    source_bytes = source.read_bytes()
    digest = hashlib.sha256(source_bytes).hexdigest()
    if digest != SOURCE_SHA256:
        raise ValueError(f"Unexpected QCElemental source SHA-256: {digest}")
    module = ast.parse(source_bytes.decode("utf-8"))
    for statement in module.body:
        if (
            isinstance(statement, ast.Assign)
            and any(
                isinstance(target, ast.Name)
                and target.id == "alvarez_2008_covalent_radii"
                for target in statement.targets
            )
        ):
            payload = ast.literal_eval(statement.value)
            if not isinstance(payload, dict):
                break
            return payload
    raise ValueError("QCElemental covalent-radius payload was not found")


def _split_label(label: str) -> tuple[str, str]:
    if "_" not in label:
        return label, "unspecified"
    symbol, variant = label.split("_", 1)
    return symbol, variant.replace("lowspin", "low_spin").replace(
        "highspin", "high_spin"
    )


def compile_source(source: Path, destination: Path) -> None:
    upstream = _source_payload(source)
    raw_records = upstream["covalent_radii"]
    records = []
    by_symbol: dict[str, list[dict[str, object]]] = {}
    for label, raw_value, comment in raw_records:
        symbol, variant = _split_label(label)
        record = {
            "symbol": symbol,
            "variant": variant,
            "value": float(raw_value),
            "comment": comment,
            "default": False,
        }
        records.append(record)
        by_symbol.setdefault(symbol, []).append(record)

    # With unknown bonding/spin, the largest published radius is the safe
    # geometric-search envelope and matches QCElemental's generic aliases.
    for variants in by_symbol.values():
        max(variants, key=lambda item: float(item["value"]))["default"] = True

    payload = {
        "metadata": {
            "dataset_id": "cristma.covalent_radii.cordero_2008",
            "version": "1",
            "original_reference": "B. Cordero et al., Dalton Trans. (2008) 2832-2838",
            "doi": "10.1039/B801115J",
            "source_repository": "https://github.com/MolSSI/QCElemental",
            "source_commit": SOURCE_COMMIT,
            "source_path": "qcelemental/data/alvarez_2008_covalent_radii.py",
            "source_license": "BSD-3-Clause",
            "source_sha256": SOURCE_SHA256,
            "generated": date.today().isoformat(),
            "record_count": len(records),
            "element_count": len(by_symbol),
            "generic_policy": "largest published variant",
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
    compile_source(arguments.source, arguments.destination)


if __name__ == "__main__":
    main()
