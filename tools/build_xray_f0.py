"""Compile the pinned xraylib FF.dat table into a CrIStMa resource."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


XRAYLIB_VERSION = "4.3.0"
XRAYLIB_COMMIT = "f94a3f5008dfd1c882b88ff26cd5052559423c83"
SOURCE_SHA256 = "9aca1801042adee51aac62ab32c9d9445e37ce5c947a7e685b42311f520c530a"
SUPPORTED_ATOMIC_NUMBERS = (1, 98)


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")


def compile_xray_f0(source: Path, destination: Path) -> None:
    raw = source.read_bytes()
    source_sha256 = hashlib.sha256(raw).hexdigest()
    if source_sha256 != SOURCE_SHA256:
        raise ValueError(
            f"unexpected FF.dat SHA-256: {source_sha256}; expected {SOURCE_SHA256}"
        )

    tokens = raw.decode("ascii").split()
    cursor = 0
    elements: list[dict[str, object]] = []
    for atomic_number in range(
        SUPPORTED_ATOMIC_NUMBERS[0],
        SUPPORTED_ATOMIC_NUMBERS[1] + 1,
    ):
        count = int(tokens[cursor])
        cursor += 1
        rows: list[tuple[float, float, float]] = []
        for _ in range(count):
            rows.append(tuple(float(item) for item in tokens[cursor : cursor + 3]))
            cursor += 3
        elements.append(
            {
                "atomic_number": atomic_number,
                "s": [row[0] for row in rows],
                "f0": [row[1] for row in rows],
                "second_derivatives": [row[2] for row in rows],
            }
        )
    if cursor != len(tokens):
        raise ValueError("FF.dat contains unexpected trailing data")

    data_sha256 = hashlib.sha256(_canonical_bytes(elements)).hexdigest()
    payload = {
        "metadata": {
            "dataset_id": "cristma.xray_f0.xraylib_ff_rayl",
            "version": "1",
            "source": "xraylib data/FF.dat",
            "xraylib_version": XRAYLIB_VERSION,
            "xraylib_commit": XRAYLIB_COMMIT,
            "source_sha256": source_sha256,
            "data_sha256": data_sha256,
            "q_convention": "sin(theta)/wavelength",
            "cristma_variable": "s = 1/(2d)",
            "s_unit": "angstrom^-1",
            "interpolation": "cubic_spline_tabulated_second_derivatives",
            "supported_atomic_numbers": list(SUPPORTED_ATOMIC_NUMBERS),
            "generator": "cristma.build_xray_f0:1",
        },
        "elements": elements,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(_canonical_bytes(payload) + b"\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    compile_xray_f0(args.source, args.destination)


if __name__ == "__main__":
    main()
