"""Compile pinned spglib tables into CRiStMa's runtime reference schema."""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
from typing import Iterable

import spglib
import spglib.error


_REQUIRED_SPGLIB_VERSION = "2.7.0"
_SCHEMA_VERSION = "1.0.0"
_DATASET_ID = "cristma.crystallography.spglib"
_TERM = re.compile(r"[+-](?:(?:\d+(?:/\d+)?)?[xyz]|\d+(?:/\d+)?)")
_CENTERING_24 = {
    "P": ((0, 0, 0),),
    "A": ((0, 0, 0), (0, 12, 12)),
    "B": ((0, 0, 0), (12, 0, 12)),
    "C": ((0, 0, 0), (12, 12, 0)),
    "I": ((0, 0, 0), (12, 12, 12)),
    "F": ((0, 0, 0), (0, 12, 12), (12, 0, 12), (12, 12, 0)),
    "R": ((0, 0, 0),),
    "H": ((0, 0, 0), (16, 8, 8), (8, 16, 16)),
}
_HEXAGONAL_R_SETTINGS = {433, 436, 444, 450, 452, 458, 460}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fraction_pair(value: Fraction) -> list[int]:
    return [value.numerator, value.denominator]


def _fraction_from_float(value: float) -> Fraction:
    result = Fraction(float(value)).limit_denominator(24)
    if abs(float(result) - float(value)) > 1e-12:
        raise ValueError(f"translation {value!r} is not representable with denominator 24")
    return result


def _crystal_system(number: int) -> str:
    if number <= 2:
        return "triclinic"
    if number <= 15:
        return "monoclinic"
    if number <= 74:
        return "orthorhombic"
    if number <= 142:
        return "tetragonal"
    if number <= 167:
        return "trigonal"
    if number <= 194:
        return "hexagonal"
    return "cubic"


def _parse_spg(path: Path) -> dict[int, list[str]]:
    records: dict[int, list[str]] = {}
    with path.open(encoding="utf-8", newline="") as source:
        for row in csv.reader(source):
            if not row:
                continue
            hall_number = int(row[0])
            if hall_number in records:
                raise ValueError(f"duplicate Hall number {hall_number}")
            if len(row) < 11:
                raise ValueError(f"incomplete spg.csv row for Hall number {hall_number}")
            records[hall_number] = row
    if not records:
        raise ValueError("spg.csv contains no records")
    return records


def _parse_component(expression: str) -> tuple[list[int], Fraction]:
    compact = expression.strip().replace(" ", "")
    signed = compact if compact.startswith(("+", "-")) else "+" + compact
    terms = _TERM.findall(signed)
    if not terms or "".join(terms) != signed:
        raise ValueError(f"unsupported Wyckoff expression: {expression!r}")
    coefficients = {"x": 0, "y": 0, "z": 0}
    translation = Fraction(0)
    for term in terms:
        sign = -1 if term[0] == "-" else 1
        body = term[1:]
        if body[-1:] in coefficients:
            variable = body[-1]
            raw_coefficient = body[:-1]
            coefficient = 1 if not raw_coefficient else int(raw_coefficient)
            coefficients[variable] += sign * coefficient
        else:
            translation += sign * Fraction(body)
    return [coefficients[axis] for axis in "xyz"], translation


def _parse_coordinate_map(source: str) -> tuple[list[list[int]], list[Fraction]]:
    text = source.strip()
    if text.startswith("(") and text.endswith(")"):
        text = text[1:-1]
    components = text.split(",")
    if len(components) != 3:
        raise ValueError(f"invalid Wyckoff coordinate: {source!r}")
    parsed = [_parse_component(component) for component in components]
    return [row for row, _translation in parsed], [translation for _row, translation in parsed]


def _position_record(
    hall_number: int,
    centering: str,
    multiplicity: int,
    letter: str,
    site_symmetry: str,
    sources: Iterable[str],
) -> dict[str, object]:
    representatives: list[dict[str, object]] = []
    shifts = _CENTERING_24[centering]
    for shift in shifts:
        for source in sources:
            matrix, translation = _parse_coordinate_map(source)
            shifted = [value + Fraction(offset, 24) for value, offset in zip(translation, shift)]
            representatives.append(
                {
                    "parameter_matrix": matrix,
                    "source": source,
                    "translation": [_fraction_pair(value % 1) for value in shifted],
                }
            )
    if len(representatives) != multiplicity:
        raise ValueError(
            f"Hall {hall_number} Wyckoff {letter} reports multiplicity {multiplicity}, "
            f"but has {len(representatives)} representatives"
        )
    return {
        "letter": letter,
        "multiplicity": multiplicity,
        "representatives": representatives,
        "site_symmetry": site_symmetry,
    }


def _parse_wyckoff(path: Path, allowed: set[int]) -> dict[str, list[dict[str, object]]]:
    groups: dict[str, list[dict[str, object]]] = {}
    current_hall: int | None = None
    current_centering: str | None = None
    current_position: dict[str, object] | None = None
    position_sources: list[str] = []

    def finish_position() -> None:
        nonlocal current_position, position_sources
        if current_position is None or current_hall is None or current_centering is None:
            return
        groups[str(current_hall)].append(
            _position_record(
                current_hall,
                current_centering,
                int(current_position["multiplicity"]),
                str(current_position["letter"]),
                str(current_position["site_symmetry"]),
                position_sources,
            )
        )
        current_position = None
        position_sources = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line == "end of data":
            break
        fields = line.split(":")
        if fields[0].isdigit():
            finish_position()
            current_hall = int(fields[0])
            if current_hall not in allowed:
                current_hall = None
                current_centering = None
                continue
            symbol = fields[1].strip()
            current_centering = "H" if current_hall in _HEXAGONAL_R_SETTINGS else symbol[0]
            if current_centering not in _CENTERING_24:
                raise ValueError(f"unsupported centering {current_centering!r}")
            groups[str(current_hall)] = []
            continue
        if current_hall is None:
            continue
        if len(fields) < 6:
            raise ValueError(f"invalid Wyckoff row: {line!r}")
        if fields[2].strip().isdigit():
            finish_position()
            current_position = {
                "multiplicity": int(fields[2]),
                "letter": fields[3].strip(),
                "site_symmetry": fields[4].strip(),
            }
        if current_position is None:
            raise ValueError(f"Wyckoff continuation without position: {line!r}")
        position_sources.extend(item.strip() for item in fields[5:9] if item.strip())
    finish_position()
    missing = allowed - {int(value) for value in groups}
    if missing:
        raise ValueError(f"Wyckoff.csv misses Hall numbers: {sorted(missing)}")
    return groups


def _metadata(
    spg_path: Path,
    wyckoff_path: Path,
    upstream_commit: str,
    compiled_date: str,
) -> dict[str, str]:
    return {
        "compiled_date": compiled_date,
        "dataset_id": _DATASET_ID,
        "license": "BSD-3-Clause",
        "schema_version": _SCHEMA_VERSION,
        "spg_sha256": _sha256(spg_path),
        "upstream": "spglib",
        "upstream_commit": upstream_commit,
        "upstream_version": _REQUIRED_SPGLIB_VERSION,
        "wyckoff_sha256": _sha256(wyckoff_path),
    }


def _write_json(path: Path, value: object) -> None:
    rendered = json.dumps(value, ensure_ascii=False, separators=(",", ":"), sort_keys=True) + "\n"
    path.write_text(rendered, encoding="utf-8")


def compile_catalog(
    spg_path: Path,
    wyckoff_path: Path,
    output_dir: Path,
    *,
    upstream_commit: str,
    compiled_date: str,
) -> tuple[Path, Path]:
    """Compile exact runtime resources from pinned spglib database inputs."""

    if spglib.__version__ != _REQUIRED_SPGLIB_VERSION:
        raise RuntimeError(f"spglib {_REQUIRED_SPGLIB_VERSION} is required")
    spglib.error.OLD_ERROR_HANDLING = False
    spg_records = _parse_spg(Path(spg_path))
    hall_numbers = set(spg_records)
    if hall_numbers != set(range(1, 531)) and not hall_numbers <= set(range(1, 531)):
        raise ValueError("Hall numbers must be a subset of 1..530")
    wyckoff_records = _parse_wyckoff(Path(wyckoff_path), hall_numbers)
    metadata = _metadata(Path(spg_path), Path(wyckoff_path), upstream_commit, compiled_date)

    groups = []
    for hall_number in sorted(hall_numbers):
        row = spg_records[hall_number]
        group_type = spglib.get_spacegroup_type(hall_number)
        symmetry = spglib.get_symmetry_from_database(hall_number)
        if group_type is None or symmetry is None:
            raise ValueError(f"spglib has no database entry for Hall number {hall_number}")
        operations = []
        for rotation, translation in zip(
            symmetry["rotations"],
            symmetry["translations"],
            strict=True,
        ):
            operations.append(
                {
                    "rotation": [[int(value) for value in axis] for axis in rotation],
                    "translation": [
                        _fraction_pair(_fraction_from_float(float(value)) % 1)
                        for value in translation
                    ],
                }
            )
        groups.append(
            {
                "centering": row[10].strip(),
                "choice": str(group_type.choice),
                "crystal_system": _crystal_system(int(group_type.number)),
                "hall_number": hall_number,
                "hall_symbol": str(group_type.hall_symbol),
                "hm_full": str(group_type.international_full),
                "hm_short": row[7].split("=")[0].strip(),
                "number": int(group_type.number),
                "operations": operations,
                "point_group": str(group_type.pointgroup_international),
            }
        )

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    groups_path = output / "space_groups.json"
    wyckoffs_path = output / "wyckoff_positions.json"
    _write_json(groups_path, {"metadata": metadata, "records": groups})
    _write_json(wyckoffs_path, {"metadata": metadata, "records": wyckoff_records})
    return groups_path, wyckoffs_path


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spg", type=Path, required=True)
    parser.add_argument("--wyckoff", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--upstream-commit", required=True)
    parser.add_argument("--compiled-date", required=True)
    args = parser.parse_args()
    compile_catalog(
        args.spg,
        args.wyckoff,
        args.output,
        upstream_commit=args.upstream_commit,
        compiled_date=args.compiled_date,
    )


if __name__ == "__main__":
    _main()
