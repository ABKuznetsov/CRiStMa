"""Complete cutoff-bounded pair search from asymmetric-unit mappings."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math

import numpy as np

from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import CrystalStructure

from .asu_mapping import AsymmetricUnitMapping, SiteImage
from .periodic_relation import PeriodicSymmetryRelation, identity_relation
from .symmetry_context import SymmetryContext, _cell_fingerprint


@dataclass(frozen=True, slots=True)
class SymmetryPairSearchPolicy:
    cutoff: float
    distance_tolerance: float = 1e-12
    max_candidates: int | None = None

    def __post_init__(self) -> None:
        for value, name, allow_zero in (
            (self.cutoff, "cutoff", False),
            (self.distance_tolerance, "distance_tolerance", True),
        ):
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise TypeError(f"{name} must be a real number")
            converted = float(value)
            if not math.isfinite(converted) or converted < 0 or (
                converted == 0 and not allow_zero
            ):
                adjective = "non-negative" if allow_zero else "positive"
                raise ValueError(f"{name} must be {adjective} and finite")
            object.__setattr__(self, name, converted)
        if self.max_candidates is not None:
            if isinstance(self.max_candidates, bool) or not isinstance(
                self.max_candidates, int
            ):
                raise TypeError("max_candidates must be an integer or None")
            if self.max_candidates <= 0:
                raise ValueError("max_candidates must be positive")


@dataclass(frozen=True, slots=True)
class SymmetryPairCandidate:
    first_site_id: str
    second_site_id: str
    relation: PeriodicSymmetryRelation
    distance: float
    vector_cartesian: tuple[float, float, float]

    def __post_init__(self) -> None:
        if not self.first_site_id or not self.second_site_id:
            raise ValueError("pair endpoint IDs must not be empty")
        if self.first_site_id > self.second_site_id:
            raise ValueError("pair endpoint IDs must use canonical order")
        if not math.isfinite(self.distance) or self.distance <= 0:
            raise ValueError("pair distance must be finite and positive")
        if len(self.vector_cartesian) != 3 or not all(
            math.isfinite(value) for value in self.vector_cartesian
        ):
            raise ValueError("pair vector must contain three finite values")


def _candidate_sort_key(candidate: SymmetryPairCandidate) -> tuple[object, ...]:
    return (
        candidate.first_site_id,
        candidate.second_site_id,
        candidate.relation,
        candidate.distance,
        candidate.vector_cartesian,
    )


@dataclass(frozen=True, slots=True)
class PairCandidateResult:
    candidates: tuple[SymmetryPairCandidate, ...]
    complete: bool
    integer_points_tested: int
    buffered_images: int
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]

    def __post_init__(self) -> None:
        if self.integer_points_tested < 0 or self.buffered_images < 0:
            raise ValueError("pair search counts must be non-negative")
        if tuple(sorted(self.candidates, key=_candidate_sort_key)) != self.candidates:
            raise ValueError("pair candidates must use deterministic order")


@dataclass(frozen=True, slots=True)
class _BufferedSiteImage:
    independent_site_id: str
    relation: PeriodicSymmetryRelation
    cartesian_position: tuple[float, float, float]


def _translation_ranges(matrix: np.ndarray, radius: float) -> tuple[range, range, range]:
    # If |v_cart| <= r and v_frac = v_cart B^-1, then
    # |v_frac_i| <= r ||column_i(B^-1)||.  The extra fractional unit covers
    # the difference of two reference-cell positions.
    component_bounds = radius * np.linalg.norm(np.linalg.inv(matrix), axis=0)
    return tuple(
        range(math.ceil(-float(bound) - 1.0), math.floor(float(bound) + 1.0) + 1)
        for bound in component_bounds
    )


def _translated_relation(
    image: SiteImage,
    translation: tuple[int, int, int],
) -> PeriodicSymmetryRelation:
    return PeriodicSymmetryRelation(
        image.representative_relation.operation_key,
        tuple(
            image.representative_relation.lattice_translation[index]
            + translation[index]
            for index in range(3)
        ),
    )


def _source_fractional(
    mapping: AsymmetricUnitMapping,
    site_id: str,
    context: SymmetryContext,
) -> tuple[float, float, float]:
    identity = identity_relation(context)
    return next(
        image.fractional_position
        for image in mapping.by_site_id[site_id].reference_cell_images
        if identity in image.equivalent_relations
    )


class SymmetryPairFinder:
    def __init__(
        self,
        cutoff: float | None = None,
        *,
        distance_tolerance: float = 1e-12,
        max_candidates: int | None = None,
        policy: SymmetryPairSearchPolicy | None = None,
    ) -> None:
        if policy is not None:
            if cutoff is not None or distance_tolerance != 1e-12 or max_candidates is not None:
                raise ValueError("provide either policy or individual search parameters")
            if not isinstance(policy, SymmetryPairSearchPolicy):
                raise TypeError("policy must be SymmetryPairSearchPolicy")
            self.policy = policy
        else:
            if cutoff is None:
                raise TypeError("cutoff is required when policy is not supplied")
            self.policy = SymmetryPairSearchPolicy(
                cutoff, distance_tolerance, max_candidates
            )

    def _validate_inputs(
        self,
        structure: CrystalStructure,
        context: SymmetryContext,
        mapping: AsymmetricUnitMapping,
    ) -> None:
        if not isinstance(structure, CrystalStructure):
            raise TypeError("structure must be CrystalStructure")
        if not isinstance(context, SymmetryContext):
            raise TypeError("context must be SymmetryContext")
        if not isinstance(mapping, AsymmetricUnitMapping):
            raise TypeError("mapping must be AsymmetricUnitMapping")
        if mapping.symmetry_context_fingerprint != context.fingerprint:
            raise ValueError("mapping belongs to another symmetry context")
        if _cell_fingerprint(structure.cell) != context.cell_fingerprint:
            raise ValueError("structure and context use different unit cells")
        if set(mapping.by_site_id) != {site.id for site in structure.sites}:
            raise ValueError("mapping and structure contain different sites")
        if not all(structure.periodic):
            raise ValueError("symmetry pair search requires a 3D periodic cell")

    def _build_buffer(
        self,
        structure: CrystalStructure,
        mapping: AsymmetricUnitMapping,
    ) -> tuple[tuple[_BufferedSiteImage, ...], int]:
        radius = self.policy.cutoff + self.policy.distance_tolerance
        ranges = _translation_ranges(structure.cell.matrix, radius)
        rows: list[_BufferedSiteImage] = []
        tested = 0
        for site_orbit in mapping.site_orbits:
            for image in site_orbit.reference_cell_images:
                for translation in product(*ranges):
                    tested += 1
                    fractional = np.asarray(image.fractional_position) + np.asarray(
                        translation, dtype=float
                    )
                    cartesian = fractional @ structure.cell.matrix
                    rows.append(
                        _BufferedSiteImage(
                            site_orbit.independent_site_id,
                            _translated_relation(image, translation),
                            tuple(float(value) for value in cartesian),
                        )
                    )
        return tuple(rows), tested

    def find_candidates(
        self,
        structure: CrystalStructure,
        context: SymmetryContext,
        mapping: AsymmetricUnitMapping,
    ) -> PairCandidateResult:
        self._validate_inputs(structure, context, mapping)
        buffer, integer_points_tested = self._build_buffer(structure, mapping)
        radius = self.policy.cutoff + self.policy.distance_tolerance
        bins: dict[tuple[int, int, int], list[_BufferedSiteImage]] = {}
        for image in buffer:
            key = tuple(math.floor(value / radius) for value in image.cartesian_position)
            bins.setdefault(key, []).append(image)

        candidates: list[SymmetryPairCandidate] = []
        tested = 0
        limit_reached = False
        coincident = False
        matrix = structure.cell.matrix
        for site_orbit in mapping.site_orbits:
            first_id = site_orbit.independent_site_id
            source = np.asarray(_source_fractional(mapping, first_id, context)) @ matrix
            source_bin = tuple(math.floor(float(value) / radius) for value in source)
            for offset in product((-1, 0, 1), repeat=3):
                key = tuple(source_bin[index] + offset[index] for index in range(3))
                for target in bins.get(key, ()):
                    if first_id > target.independent_site_id:
                        continue
                    if self.policy.max_candidates is not None and tested >= self.policy.max_candidates:
                        limit_reached = True
                        break
                    tested += 1
                    vector = np.asarray(target.cartesian_position) - source
                    distance = float(np.linalg.norm(vector))
                    if distance <= self.policy.distance_tolerance:
                        is_own_stabilizer = (
                            first_id == target.independent_site_id
                            and target.relation
                            in mapping.by_site_id[first_id].stabilizer_relations
                        )
                        coincident = coincident or not is_own_stabilizer
                    elif distance <= radius:
                        candidates.append(
                            SymmetryPairCandidate(
                                first_id,
                                target.independent_site_id,
                                target.relation,
                                distance,
                                tuple(float(value) for value in vector),
                            )
                        )
                if limit_reached:
                    break
            if limit_reached:
                break

        diagnostics: list[Diagnostic] = []
        if limit_reached:
            diagnostics.append(
                Diagnostic(
                    Severity.WARNING,
                    "symmetry.pairs.search_limit_reached",
                    "pair search stopped at the configured candidate limit",
                    recovery="increase max_candidates to cover the complete cutoff region",
                )
            )
        if coincident:
            diagnostics.append(
                Diagnostic(
                    Severity.WARNING,
                    "symmetry.pairs.coincident_sites_excluded",
                    "zero-distance site images were excluded from geometric contacts",
                )
            )
        return PairCandidateResult(
            tuple(sorted(candidates, key=_candidate_sort_key)),
            not limit_reached,
            integer_points_tested,
            len(buffer),
            tuple(diagnostics),
            (
                ("method", "asu_cartesian_spatial_bins"),
                ("cutoff", self.policy.cutoff),
                ("distance_tolerance", self.policy.distance_tolerance),
                ("max_candidates", self.policy.max_candidates),
                ("candidate_pairs_tested", tested),
            ),
        )


__all__ = [
    "PairCandidateResult",
    "SymmetryPairCandidate",
    "SymmetryPairFinder",
    "SymmetryPairSearchPolicy",
]
