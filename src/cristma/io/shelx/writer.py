"""Preserving source and canonical scientific SHELX writers."""

from __future__ import annotations

from dataclasses import dataclass
import math

from cristma.core.values import MeasuredValue
from cristma.structure import CrystalStructure, IndependentSite
from cristma.symmetry.affine import AffineOperation

from .document import ShelxDocument
from .occupancy import ShelxOccupancyExpression
from .symmetry import format_shelx_symmetry


@dataclass(frozen=True, slots=True)
class ShelxWriteOptions:
    """Measurement information and output choices required by SHELX."""

    wavelength: float | MeasuredValue | None = None
    hklf: int = 4
    title: str | None = None

    def __post_init__(self) -> None:
        numeric = (
            self.wavelength.value
            if isinstance(self.wavelength, MeasuredValue)
            else self.wavelength
        )
        if numeric is not None and (not math.isfinite(numeric) or numeric <= 0):
            raise ValueError("wavelength must be a positive finite number")
        if isinstance(self.hklf, bool) or not isinstance(self.hklf, int):
            raise TypeError("hklf must be an integer")


def _number(value: MeasuredValue | float | int) -> str:
    if isinstance(value, MeasuredValue):
        if value.raw not in {None, ""}:
            return str(value.raw)
        if value.value is None:
            raise ValueError("canonical SHELX output cannot write a missing number")
        value = value.value
    return format(float(value), ".15g")


def _is_identity(operation: AffineOperation) -> bool:
    return format_shelx_symmetry(operation.normalized()) == "x,y,z"


def _site_labels(site: IndependentSite) -> tuple[str, ...]:
    if len(site.components) == 1:
        return (site.label,)
    return tuple(
        f"{site.label}_{component.element or index}"
        for index, component in enumerate(site.components, start=1)
    )


def _occupancy(component: object) -> str:
    metadata = getattr(component, "metadata")
    expression = metadata.get("shelx_occupancy")
    if isinstance(expression, ShelxOccupancyExpression):
        return expression.raw
    return _number(10.0 + float(getattr(component, "occupancy").value))


def _displacement(site: IndependentSite) -> tuple[str, ...]:
    displacement = site.displacement
    if displacement is None:
        return ("0.05",)
    if displacement.kind in {"U_iso", "B_iso"}:
        if displacement.isotropic is None:
            raise ValueError(f"site {site.label!r} has no isotropic displacement value")
        return (_number(displacement.isotropic),)
    if displacement.kind == "U_aniso" and displacement.tensor is not None:
        tensor = displacement.tensor
        return tuple(
            _number(value)
            for value in (
                tensor[0][0],
                tensor[1][1],
                tensor[2][2],
                tensor[1][2],
                tensor[0][2],
                tensor[0][1],
            )
        )
    raise ValueError(f"unsupported displacement model for site {site.label!r}")


def _elements(crystal: CrystalStructure) -> tuple[str, ...]:
    result: list[str] = []
    for site in crystal.sites:
        for component in site.components:
            element = component.element
            if element is None:
                raise ValueError("canonical SHELX output requires known elements")
            if element not in result:
                result.append(element)
    return tuple(result)


def _unit_contents(
    crystal: CrystalStructure,
    elements: tuple[str, ...],
) -> tuple[float, ...]:
    totals = {element: 0.0 for element in elements}
    for atom in crystal.atomic_view().atoms:
        for component in atom.components:
            element = component.element
            if element in totals:
                totals[element] += float(component.occupancy.value)
    return tuple(totals[element] for element in elements)


def write_shelx_document(
    document: ShelxDocument,
    *,
    mode: str = "preserve",
) -> str:
    """Render a SHELX document without altering untouched source text."""

    if mode != "preserve":
        raise ValueError("ShelxDocument supports only preserve-mode writing")
    return document.render_preserved()


def write_crystal_shelx(
    crystal: CrystalStructure,
    *,
    options: ShelxWriteOptions | None = None,
) -> str:
    """Write a canonical SHELX instruction file from a crystal snapshot."""

    if options is None or options.wavelength is None:
        raise ValueError("canonical SHELX writing requires a wavelength")
    wavelength = options.wavelength
    if isinstance(wavelength, MeasuredValue) and wavelength.value is None:
        raise ValueError("canonical SHELX writing requires a wavelength")

    title = options.title or crystal.metadata.get("shelx_title") or crystal.name
    cell = crystal.cell
    lines = [
        f"TITL {title}",
        "CELL " + " ".join(
            _number(value)
            for value in (
                wavelength,
                cell.a,
                cell.b,
                cell.c,
                cell.alpha,
                cell.beta,
                cell.gamma,
            )
        ),
        "LATT -1",
    ]
    if crystal.space_group is not None:
        lines.extend(
            f"SYMM {format_shelx_symmetry(operation)}"
            for operation in crystal.space_group.operations
            if not _is_identity(operation)
        )

    elements = _elements(crystal)
    lines.append("SFAC " + " ".join(elements))
    lines.append(
        "UNIT " + " ".join(_number(value) for value in _unit_contents(crystal, elements))
    )
    free_variables = tuple(crystal.metadata.get("shelx_free_variables", ()))
    if free_variables:
        lines.append("FVAR " + " ".join(_number(value) for value in free_variables))

    sfac_indices = {element: index for index, element in enumerate(elements, start=1)}
    for site in crystal.sites:
        coordinates = " ".join(_number(value) for value in site.fractional)
        displacement = " ".join(_displacement(site))
        for label, component in zip(_site_labels(site), site.components, strict=True):
            lines.append(
                f"{label} {sfac_indices[component.element]} {coordinates} "
                f"{_occupancy(component)} {displacement}"
            )
    lines.extend((f"HKLF {options.hklf}", "END"))
    return "\n".join(lines) + "\n"


__all__ = ["ShelxWriteOptions", "write_crystal_shelx", "write_shelx_document"]
