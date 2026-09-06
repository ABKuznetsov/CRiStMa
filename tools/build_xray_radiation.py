"""Compile pinned xraylib Cu K-alpha lines into a CrIStMa resource."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


XRAYLIB_VERSION = "4.3.0"
XRAYLIB_COMMIT = "f94a3f5008dfd1c882b88ff26cd5052559423c83"
HC_KEV_ANGSTROM = 12.398419843320026
COPPER_ATOMIC_NUMBER = 29


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")


def compile_xray_radiation(xraylib: Any, destination: Path) -> None:
    """Extract Cu K-alpha line energies and radiative rates."""

    reported_version = str(xraylib.__version__)
    if reported_version != XRAYLIB_VERSION:
        raise ValueError(
            f"unexpected xraylib version: {reported_version}; expected {XRAYLIB_VERSION}"
        )
    components: list[dict[str, object]] = []
    for component_id, label, line in (
        ("cu-ka1", "Cu Kalpha1", xraylib.KA1_LINE),
        ("cu-ka2", "Cu Kalpha2", xraylib.KA2_LINE),
    ):
        energy_kev = float(xraylib.LineEnergy(COPPER_ATOMIC_NUMBER, line))
        relative_weight = float(xraylib.RadRate(COPPER_ATOMIC_NUMBER, line))
        if energy_kev <= 0.0 or relative_weight <= 0.0:
            raise ValueError(f"xraylib returned invalid data for {component_id}")
        components.append(
            {
                "component_id": component_id,
                "label": label,
                "energy_kev": energy_kev,
                "wavelength_angstrom": HC_KEV_ANGSTROM / energy_kev,
                "relative_weight": relative_weight,
            }
        )
    resource_checksum = hashlib.sha256(_canonical_bytes(components)).hexdigest()
    payload = {
        "metadata": {
            "dataset_id": "cristma.xray_radiation.xraylib_cu_ka",
            "version": "1",
            "source": "xraylib Cu line data",
            "energy_source": "xraylib LineEnergy",
            "radiative_rate_source": "xraylib RadRate",
            "energy_to_wavelength_formula": (
                "lambda_angstrom = hc_keV_angstrom / energy_keV"
            ),
            "hc_value": HC_KEV_ANGSTROM,
            "hc_units": "keV angstrom",
            "xraylib_version": XRAYLIB_VERSION,
            "xraylib_commit": XRAYLIB_COMMIT,
            "resource_checksum": resource_checksum,
            "generator": "cristma.build_xray_radiation:1",
        },
        "components": components,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(_canonical_bytes(payload) + b"\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    import xraylib

    compile_xray_radiation(xraylib, args.destination)


if __name__ == "__main__":
    main()
