"""Compile pinned xraylib laboratory K-alpha lines into a CrIStMa resource."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


XRAYLIB_VERSION = "4.3.0"
XRAYLIB_COMMIT = "f94a3f5008dfd1c882b88ff26cd5052559423c83"
HC_KEV_ANGSTROM = 12.398419843320026
TUBE_TARGETS = (
    ("Cr", 24),
    ("Fe", 26),
    ("Co", 27),
    ("Cu", 29),
    ("Mo", 42),
    ("Ag", 47),
)


def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("ascii")


def compile_xray_radiation(xraylib: Any, destination: Path) -> None:
    """Extract common tube-target K-alpha energies and radiative rates."""

    reported_version = str(xraylib.__version__)
    if reported_version != XRAYLIB_VERSION:
        raise ValueError(
            f"unexpected xraylib version: {reported_version}; expected {XRAYLIB_VERSION}"
        )
    spectra: list[dict[str, object]] = []
    for symbol, atomic_number in TUBE_TARGETS:
        components: list[dict[str, object]] = []
        for suffix, line in (
            ("ka1", xraylib.KA1_LINE),
            ("ka2", xraylib.KA2_LINE),
        ):
            component_id = f"{symbol.lower()}-{suffix}"
            energy_kev = float(xraylib.LineEnergy(atomic_number, line))
            relative_weight = float(xraylib.RadRate(atomic_number, line))
            if energy_kev <= 0.0 or relative_weight <= 0.0:
                raise ValueError(f"xraylib returned invalid data for {component_id}")
            components.append(
                {
                    "component_id": component_id,
                    "label": f"{symbol} Kalpha{suffix[-1]}",
                    "energy_kev": energy_kev,
                    "wavelength_angstrom": HC_KEV_ANGSTROM / energy_kev,
                    "relative_weight": relative_weight,
                }
            )
        spectra.append(
            {
                "source_id": f"xray-tube:{symbol.lower()}-ka",
                "target_element": symbol,
                "atomic_number": atomic_number,
                "components": components,
            }
        )
    resource_checksum = hashlib.sha256(_canonical_bytes(spectra)).hexdigest()
    payload = {
        "metadata": {
            "dataset_id": "cristma.xray_radiation.xraylib_tube_k_alpha",
            "version": "2",
            "source": "xraylib laboratory tube K-alpha line data",
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
            "generator": "cristma.build_xray_radiation:2",
        },
        "spectra": spectra,
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
