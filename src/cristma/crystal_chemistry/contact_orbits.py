"""Symmetry-equivalent families of resolved periodic contacts."""

from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
import math

from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import AtomicView, ExpandedAtom
from cristma.symmetry import AffineOperation

from .contacts import ResolutionStatus, ResolvedContact, ResolvedContactOrbit


SYMMETRY_INCOMPLETE = "crystal_chemistry.contact_orbit.symmetry_incomplete"
CLASSIFICATION_INCONSISTENT = (
    "crystal_chemistry.contact_orbit.classification_inconsistent"
)

Translation = tuple[int, int, int]
ContactKey = tuple[str, str, Translation]
MemberDescriptor = tuple[str, str, str, str, Translation]


@dataclass(frozen=True, slots=True)
class ContactOrbitBuildResult:
    contacts: tuple[ResolvedContact, ...]
    contact_orbits: tuple[ResolvedContactOrbit, ...]
    diagnostics: tuple[Diagnostic, ...]
    complete: bool


def aggregate_resolution_status(
    statuses: tuple[ResolutionStatus, ...],
    *,
    applicable: bool,
) -> ResolutionStatus:
    """Aggregate scientific outcomes without interpreting diagnostic prose."""

    if ResolutionStatus.INCOMPLETE in statuses:
        return ResolutionStatus.INCOMPLETE
    if ResolutionStatus.AMBIGUOUS in statuses:
        return ResolutionStatus.AMBIGUOUS
    return ResolutionStatus.RESOLVED if applicable else ResolutionStatus.NOT_APPLICABLE


def _negative(value: Translation) -> Translation:
    return (-value[0], -value[1], -value[2])


def _canonical_key(first: str, second: str, translation: Translation) -> ContactKey:
    forward = (first, second, translation)
    reverse = (second, first, _negative(translation))
    return min(forward, reverse)


def _wrap(value: tuple[float, float, float], tolerance: float):
    wrapped: list[float] = []
    lattice: list[int] = []
    for item in value:
        nearest = round(item)
        normalized = float(nearest) if abs(item - nearest) <= tolerance else item
        shift = math.floor(normalized)
        coordinate = normalized - shift
        if abs(coordinate - 1.0) <= tolerance:
            coordinate = 0.0
            shift += 1
        wrapped.append(0.0 if abs(coordinate) <= tolerance else coordinate)
        lattice.append(int(shift))
    return tuple(wrapped), tuple(lattice)


def _transform(
    operation: AffineOperation,
    fractional: tuple[float, float, float],
    cell_translation: Translation,
    tolerance: float,
):
    global_position = tuple(
        coordinate + translation
        for coordinate, translation in zip(
            fractional, cell_translation, strict=True
        )
    )
    raw = tuple(
        math.fsum(
            float(coefficient) * coordinate
            for coefficient, coordinate in zip(row, global_position, strict=True)
        )
        + float(offset)
        for row, offset in zip(
            operation.rotation, operation.translation, strict=True
        )
    )
    return _wrap(raw, tolerance)


def _periodically_equal(
    left: tuple[float, float, float],
    right: tuple[float, float, float],
    tolerance: float,
) -> bool:
    return all(
        abs((a - b + 0.5) % 1.0 - 0.5) <= tolerance + 1e-12
        for a, b in zip(left, right, strict=True)
    )


def _match_atom(
    atoms_by_site: dict[str, tuple[ExpandedAtom, ...]],
    source_site_id: str,
    fractional: tuple[float, float, float],
    tolerance: float,
) -> ExpandedAtom | None:
    matches = tuple(
        atom
        for atom in atoms_by_site.get(source_site_id, ())
        if _periodically_equal(atom.fractional, fractional, tolerance)
    )
    return min(matches, key=lambda item: item.id) if matches else None


def _contact_key(contact: ResolvedContact) -> ContactKey:
    geometric = contact.geometric_contact
    return _canonical_key(
        geometric.first_atom_id,
        geometric.second_atom_id,
        geometric.cell_translation or (0, 0, 0),
    )


def _image_descriptors(atom: ExpandedAtom):
    return tuple(
        sorted(
            {
                (image.operation_id, image.normalization_translation)
                for image in atom.equivalent_images
            }
        )
    )


def _member_descriptor(
    contact: ResolvedContact,
    atoms: dict[str, ExpandedAtom],
) -> MemberDescriptor:
    geometric = contact.geometric_contact
    first = atoms[geometric.first_atom_id]
    second = atoms[geometric.second_atom_id]
    translation = geometric.cell_translation or (0, 0, 0)
    candidates: list[MemberDescriptor] = []
    for first_operation, first_normalization in _image_descriptors(first):
        for second_operation, second_normalization in _image_descriptors(second):
            relative = tuple(
                second_normalization[index]
                + translation[index]
                - first_normalization[index]
                for index in range(3)
            )
            forward: MemberDescriptor = (
                first.source_site_id,
                first_operation,
                second.source_site_id,
                second_operation,
                relative,
            )
            reverse: MemberDescriptor = (
                second.source_site_id,
                second_operation,
                first.source_site_id,
                first_operation,
                _negative(relative),
            )
            candidates.append(min(forward, reverse))
    return min(candidates)


def _orbit_id(descriptors: tuple[MemberDescriptor, ...]) -> str:
    payload = json.dumps(
        descriptors,
        ensure_ascii=True,
        separators=(",", ":"),
    ).encode("ascii")
    return f"contact-orbit:{hashlib.sha256(payload).hexdigest()}"


def build_contact_orbits(
    view: AtomicView[ExpandedAtom],
    contacts: tuple[ResolvedContact, ...],
    operations: tuple[AffineOperation, ...],
    tolerance: float,
) -> ContactOrbitBuildResult:
    """Group resolved contacts by the exact space-group action on endpoints."""

    if view.fractional is None or not all(view.periodic):
        raise ValueError("contact orbit construction requires a 3D periodic atomic view")
    if not operations:
        raise ValueError("contact orbit construction requires symmetry operations")
    if not math.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("contact orbit tolerance must be positive and finite")
    atoms = {atom.id: atom for atom in view.atoms}
    atoms_by_site: dict[str, list[ExpandedAtom]] = {}
    for atom in view.atoms:
        atoms_by_site.setdefault(atom.source_site_id, []).append(atom)
    frozen_atoms_by_site = {
        key: tuple(sorted(values, key=lambda item: item.id))
        for key, values in atoms_by_site.items()
    }
    by_key: dict[ContactKey, ResolvedContact] = {}
    for contact in contacts:
        key = _contact_key(contact)
        if key in by_key:
            raise ValueError("resolved contacts must identify unique physical contacts")
        if key[0] not in atoms or key[1] not in atoms:
            raise ValueError("resolved contact references an atom outside the atomic view")
        by_key[key] = contact

    parent = {key: key for key in by_key}

    def find(key: ContactKey) -> ContactKey:
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(left: ContactKey, right: ContactKey) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    diagnostics: list[Diagnostic] = []
    complete = True
    for key, contact in by_key.items():
        geometric = contact.geometric_contact
        first = atoms[geometric.first_atom_id]
        second = atoms[geometric.second_atom_id]
        translation = geometric.cell_translation or (0, 0, 0)
        missing = False
        for operation in operations:
            first_fractional, first_lattice = _transform(
                operation, first.fractional, (0, 0, 0), tolerance
            )
            second_fractional, second_lattice = _transform(
                operation, second.fractional, translation, tolerance
            )
            first_image = _match_atom(
                frozen_atoms_by_site,
                first.source_site_id,
                first_fractional,
                tolerance,
            )
            second_image = _match_atom(
                frozen_atoms_by_site,
                second.source_site_id,
                second_fractional,
                tolerance,
            )
            if first_image is None or second_image is None:
                missing = True
                continue
            transformed_translation = tuple(
                second_lattice[index] - first_lattice[index]
                for index in range(3)
            )
            transformed_key = _canonical_key(
                first_image.id,
                second_image.id,
                transformed_translation,
            )
            if transformed_key not in by_key:
                missing = True
                continue
            union(key, transformed_key)
        if missing:
            complete = False
            diagnostics.append(
                Diagnostic(
                    Severity.WARNING,
                    SYMMETRY_INCOMPLETE,
                    f"Contact {geometric.contact_id} is missing a resolved symmetry image",
                )
            )

    groups: dict[ContactKey, list[ResolvedContact]] = {}
    for key, contact in by_key.items():
        groups.setdefault(find(key), []).append(contact)

    replacements: dict[str, ResolvedContact] = {}
    orbit_rows: list[ResolvedContactOrbit] = []
    for members in groups.values():
        described = tuple(
            sorted(
                ((_member_descriptor(item, atoms), item) for item in members),
                key=lambda row: (
                    row[0], row[1].geometric_contact.contact_id
                ),
            )
        )
        descriptors = tuple(row[0] for row in described)
        contact_orbit_id = _orbit_id(descriptors)
        replaced = tuple(
            replace(item, contact_orbit_id=contact_orbit_id)
            for _, item in described
        )
        classifications = {item.contact_classification for item in replaced}
        if len(classifications) != 1:
            complete = False
            diagnostics.append(
                Diagnostic(
                    Severity.WARNING,
                    CLASSIFICATION_INCONSISTENT,
                    "Symmetry-equivalent contacts have inconsistent classifications",
                )
            )
        for item in replaced:
            replacements[item.geometric_contact.contact_id] = item
        orbit_rows.append(
            ResolvedContactOrbit(
                contact_orbit_id,
                replaced[0].geometric_contact.contact_id,
                replaced,
            )
        )

    contact_orbits = tuple(
        sorted(orbit_rows, key=lambda item: item.contact_orbit_id)
    )
    output_contacts = tuple(
        replacements[item.geometric_contact.contact_id] for item in contacts
    )
    unique_diagnostics = tuple(dict.fromkeys(diagnostics))
    return ContactOrbitBuildResult(
        output_contacts,
        contact_orbits,
        unique_diagnostics,
        complete,
    )


__all__ = [
    "CLASSIFICATION_INCONSISTENT",
    "ContactOrbitBuildResult",
    "SYMMETRY_INCOMPLETE",
    "aggregate_resolution_status",
    "build_contact_orbits",
]
