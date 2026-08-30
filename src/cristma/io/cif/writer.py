"""Preserving and canonical CIF 1.1 writers."""

from __future__ import annotations

from fractions import Fraction
import re

from cristma.structure import CrystalStructure, IndependentSite
from cristma.core.values import MeasuredValue
from cristma.symmetry.affine import AffineOperation

from .document import CifDocument


def write_cif_document(
    document: CifDocument,
    *,
    mode: str = "preserve",
) -> str:
    """Render a parsed CIF document while retaining untouched source text."""

    if mode != "preserve":
        raise ValueError("CifDocument supports only preserve-mode writing")
    if not document.edits:
        return document.raw_source

    ordered = sorted(document.edits, key=lambda item: (item.start, item.end))
    for previous, current in zip(ordered, ordered[1:]):
        if current.start < previous.end:
            raise ValueError("overlapping CIF source edits")

    rendered = document.raw_source
    for edit in reversed(ordered):
        rendered = rendered[: edit.start] + edit.replacement + rendered[edit.end :]
    return rendered


def _number(value: MeasuredValue) -> str:
    if value.raw not in {None, ""}:
        return str(value.raw)
    if value.value is None:
        if value.missing.value == "unknown":
            return "?"
        if value.missing.value == "inapplicable":
            return "."
        return "?"
    return format(float(value.value), ".15g")


def _quote(value: str) -> str:
    if "\n" in value or "\r" in value:
        normalized = value.replace("\r\n", "\n").replace("\r", "\n")
        return f";{normalized}\n;"
    folded = value.casefold()
    reserved = (
        value == ""
        or value[0] in "_#$;"
        or any(character.isspace() for character in value)
        or "," in value
        or folded in {"loop_", "stop_", "global_"}
        or folded.startswith(("data_", "save_"))
    )
    if not reserved:
        return value
    if "'" not in value:
        return f"'{value}'"
    if '"' not in value:
        return f'"{value}"'
    return f";{value}\n;"


def _fraction(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _operation_component(
    row: tuple[Fraction, Fraction, Fraction],
    translation: Fraction,
) -> str:
    terms: list[str] = []
    for coefficient, variable in zip(row, "xyz", strict=True):
        if coefficient == 0:
            continue
        if coefficient not in {Fraction(-1), Fraction(1)}:
            raise ValueError("canonical CIF symmetry supports unit variable coefficients")
        if coefficient < 0:
            terms.append(f"-{variable}")
        else:
            terms.append(("+" if terms else "") + variable)
    normalized = translation % 1
    if normalized:
        terms.append(("+" if terms else "") + _fraction(normalized))
    return "".join(terms) or "0"


def _operation(operation: AffineOperation) -> str:
    if operation.source:
        return operation.source
    return ",".join(
        _operation_component(row, offset)
        for row, offset in zip(
            operation.rotation,
            operation.translation,
            strict=True,
        )
    )


def _block_name(value: str) -> str:
    normalized = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    return normalized or "cristma"


def _component_labels(site: IndependentSite) -> tuple[str, ...]:
    reported = tuple(site.label.split("/"))
    if len(reported) == len(site.components):
        return reported
    if len(site.components) == 1:
        return (site.label,)
    return tuple(
        f"{component.element}{index}"
        for index, component in enumerate(site.components, start=1)
    )


def _metadata_lines(crystal: CrystalStructure) -> list[str]:
    tag_by_key = {
        "mineral_name": "_chemical_name_mineral",
        "common_name": "_chemical_name_common",
        "systematic_name": "_chemical_name_systematic",
        "publication_title": "_publ_section_title",
        "journal": "_journal_name_full",
        "year": "_journal_year",
        "volume": "_journal_volume",
        "page_first": "_journal_page_first",
        "page_last": "_journal_page_last",
        "doi": "_journal_paper_doi",
    }
    return [
        f"{tag_by_key[key]} {_quote(str(crystal.metadata[key]))}"
        for key in tag_by_key
        if key in crystal.metadata
    ]


def _atom_headers(crystal: CrystalStructure) -> list[tuple[str, str]]:
    headers = [
        ("label", "_atom_site_label"),
        ("element", "_atom_site_type_symbol"),
        ("x", "_atom_site_fract_x"),
        ("y", "_atom_site_fract_y"),
        ("z", "_atom_site_fract_z"),
        ("occupancy", "_atom_site_occupancy"),
    ]
    if any(component.oxidation_state is not None for site in crystal.sites for component in site.components):
        headers.append(("oxidation", "_atom_site_oxidation_number"))
    if any(site.wyckoff is not None for site in crystal.sites):
        headers.append(("wyckoff", "_atom_site_Wyckoff_symbol"))
    if any(
        site.reported_multiplicity is not None or site.calculated_multiplicity is not None
        for site in crystal.sites
    ):
        headers.append(("multiplicity", "_atom_site_symmetry_multiplicity"))
    if any(site.disorder_assembly is not None for site in crystal.sites):
        headers.append(("assembly", "_atom_site_disorder_assembly"))
    if any(site.disorder_group is not None for site in crystal.sites):
        headers.append(("group", "_atom_site_disorder_group"))
    if any(site.displacement is not None and site.displacement.kind == "U_iso" for site in crystal.sites):
        headers.append(("u_iso", "_atom_site_U_iso_or_equiv"))
    if any(site.displacement is not None and site.displacement.kind == "B_iso" for site in crystal.sites):
        headers.append(("b_iso", "_atom_site_B_iso_or_equiv"))
    return headers


def _atom_rows(
    crystal: CrystalStructure,
    headers: list[tuple[str, str]],
) -> list[list[str]]:
    rows: list[list[str]] = []
    for site in crystal.sites:
        labels = _component_labels(site)
        for label, component in zip(labels, site.components, strict=True):
            values = {
                "label": label,
                "element": component.element,
                "x": _number(site.fractional[0]),
                "y": _number(site.fractional[1]),
                "z": _number(site.fractional[2]),
                "occupancy": _number(component.occupancy),
                "oxidation": (
                    _number(component.oxidation_state)
                    if component.oxidation_state is not None
                    else "?"
                ),
                "wyckoff": site.wyckoff or "?",
                "multiplicity": str(
                    site.reported_multiplicity
                    if site.reported_multiplicity is not None
                    else site.calculated_multiplicity
                ) if (
                    site.reported_multiplicity is not None
                    or site.calculated_multiplicity is not None
                ) else "?",
                "assembly": site.disorder_assembly or ".",
                "group": site.disorder_group or ".",
                "u_iso": (
                    _number(site.displacement.isotropic)
                    if site.displacement is not None
                    and site.displacement.kind == "U_iso"
                    and site.displacement.isotropic is not None
                    else "?"
                ),
                "b_iso": (
                    _number(site.displacement.isotropic)
                    if site.displacement is not None
                    and site.displacement.kind == "B_iso"
                    and site.displacement.isotropic is not None
                    else "?"
                ),
            }
            rows.append([_quote(values[key]) for key, _tag in headers])
    return rows


def _anisotropic_lines(crystal: CrystalStructure) -> list[str]:
    sites = [
        site
        for site in crystal.sites
        if site.displacement is not None
        and site.displacement.kind == "U_aniso"
        and site.displacement.tensor is not None
    ]
    if not sites:
        return []
    lines = [
        "loop_",
        "_atom_site_aniso_label",
        "_atom_site_aniso_U_11",
        "_atom_site_aniso_U_22",
        "_atom_site_aniso_U_33",
        "_atom_site_aniso_U_12",
        "_atom_site_aniso_U_13",
        "_atom_site_aniso_U_23",
    ]
    for site in sites:
        tensor = site.displacement.tensor
        values = (
            tensor[0][0],
            tensor[1][1],
            tensor[2][2],
            tensor[0][1],
            tensor[0][2],
            tensor[1][2],
        )
        for label in _component_labels(site):
            lines.append(" ".join([_quote(label), *(_number(value) for value in values)]))
    return lines


def write_crystal_cif(
    crystal: CrystalStructure,
    *,
    block_name: str | None = None,
) -> str:
    """Write a normalized CIF containing the canonical asymmetric unit."""

    lines = [f"data_{_block_name(block_name or crystal.name)}"]
    lines.extend(
        (
            f"_cell_length_a {_number(crystal.cell.a)}",
            f"_cell_length_b {_number(crystal.cell.b)}",
            f"_cell_length_c {_number(crystal.cell.c)}",
            f"_cell_angle_alpha {_number(crystal.cell.alpha)}",
            f"_cell_angle_beta {_number(crystal.cell.beta)}",
            f"_cell_angle_gamma {_number(crystal.cell.gamma)}",
        )
    )
    group = crystal.space_group
    if group is not None:
        if group.hm_symbol:
            lines.append(f"_space_group_name_H-M_alt {_quote(group.hm_symbol)}")
        if group.hall_symbol:
            lines.append(f"_space_group_name_Hall {_quote(group.hall_symbol)}")
        if group.number is not None:
            lines.append(f"_space_group_IT_number {group.number}")
    if crystal.formula:
        lines.append(f"_chemical_formula_sum {_quote(crystal.formula)}")
    lines.extend(_metadata_lines(crystal))

    if group is not None:
        lines.extend(("loop_", "_space_group_symop_operation_xyz"))
        lines.extend(_quote(_operation(operation)) for operation in group.operations)

    headers = _atom_headers(crystal)
    lines.append("loop_")
    lines.extend(tag for _key, tag in headers)
    lines.extend(" ".join(row) for row in _atom_rows(crystal, headers))
    lines.extend(_anisotropic_lines(crystal))
    return "\n".join(lines) + "\n"
