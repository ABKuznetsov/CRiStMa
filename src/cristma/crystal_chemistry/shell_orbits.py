"""Multiplicity-weighted coordination shells built from incidence orbits."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
import math

from cristma.crystallography import PairTableStatus, SymmetryPairTable
from cristma.crystallography.symmetry_context import _digest
from cristma.diagnostics import Diagnostic, Severity

from .contacts import EvidenceStatus, ResolutionStatus, SecondaryEvidence
from .incidence_orbits import ContactIncidenceOrbit
from .orbit_contacts import OrientationMode, ResolvedContactOrbit
from .policy import ShellResolutionPolicy


class ShellRole(StrEnum):
    PRIMARY = "primary"
    SECONDARY = "secondary"


def _interpretation_context_id(interpretation) -> str:
    return "shell-interpretation-context:" + _digest(
        (
            interpretation.interaction_type.value,
            interpretation.interaction_layer.value,
            interpretation.grammar_priority.value,
        )
    )


@dataclass(frozen=True, slots=True)
class CoordinationShellAlternative:
    alternative_id: str
    primary_incidence_ids: tuple[str, ...]
    secondary_incidence_ids: tuple[str, ...]
    geometric_CN: int
    mean_occupied_neighbors: float
    boundary_evidence: tuple[SecondaryEvidence, ...]
    status: ResolutionStatus

    def __post_init__(self) -> None:
        if not self.alternative_id:
            raise ValueError("shell alternative ID must not be empty")
        if tuple(sorted(set(self.primary_incidence_ids))) != self.primary_incidence_ids:
            raise ValueError("primary incidence IDs must be unique and sorted")
        if tuple(sorted(set(self.secondary_incidence_ids))) != self.secondary_incidence_ids:
            raise ValueError("secondary incidence IDs must be unique and sorted")
        if set(self.primary_incidence_ids) & set(self.secondary_incidence_ids):
            raise ValueError("primary and secondary incidence sets must be disjoint")
        if self.geometric_CN < 0:
            raise ValueError("geometric coordination number must be non-negative")
        if not math.isfinite(self.mean_occupied_neighbors) or self.mean_occupied_neighbors < 0:
            raise ValueError("occupied-neighbour count must be non-negative and finite")


@dataclass(frozen=True, slots=True)
class CoordinationShellOrbit:
    shell_orbit_id: str
    center_independent_site_id: str
    interpretation_context_id: str
    selected_alternative: str | None
    alternatives: tuple[CoordinationShellAlternative, ...]
    status: ResolutionStatus
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]

    def __post_init__(self) -> None:
        if not self.shell_orbit_id or not self.center_independent_site_id:
            raise ValueError("coordination-shell identities must not be empty")
        if not self.interpretation_context_id:
            raise ValueError("shell interpretation context must not be empty")
        ids = tuple(item.alternative_id for item in self.alternatives)
        if len(ids) != len(set(ids)):
            raise ValueError("shell alternative IDs must be unique")
        if self.status is ResolutionStatus.RESOLVED:
            if self.selected_alternative is None or self.selected_alternative not in ids:
                raise ValueError("resolved shell must select exactly one alternative")
        elif self.selected_alternative is not None:
            raise ValueError("unresolved shell must not select an alternative")
        if self.status is ResolutionStatus.AMBIGUOUS and len(self.alternatives) < 2:
            raise ValueError("ambiguous shell requires at least two alternatives")
        if self.status is ResolutionStatus.NOT_APPLICABLE:
            raise ValueError("non-applicable coordination requests do not create shell orbits")

    @property
    def selected(self) -> CoordinationShellAlternative | None:
        if self.selected_alternative is None:
            return None
        return next(
            item for item in self.alternatives
            if item.alternative_id == self.selected_alternative
        )


@dataclass(frozen=True, slots=True)
class _IncidenceDistance:
    incidence: ContactIncidenceOrbit
    lower: float
    central: float
    upper: float


@dataclass(frozen=True, slots=True)
class _BoundaryCandidate:
    primary_ids: tuple[str, ...]
    relative_gap: float
    internal_spread: float


@dataclass(frozen=True, slots=True)
class _BoundaryDecision:
    status: ResolutionStatus
    selected: _BoundaryCandidate | None
    alternatives: tuple[_BoundaryCandidate, ...]
    diagnostics: tuple[Diagnostic, ...]


def _distance_groups(
    records: tuple[_IncidenceDistance, ...],
    attribute: str,
    tolerance: float,
) -> tuple[tuple[_IncidenceDistance, ...], ...]:
    ordered = tuple(sorted(records, key=lambda item: (getattr(item, attribute), item.incidence.incidence_orbit_id)))
    groups: list[list[_IncidenceDistance]] = []
    for record in ordered:
        if not groups or getattr(record, attribute) - getattr(groups[-1][-1], attribute) > tolerance + 1e-12:
            groups.append([record])
        else:
            groups[-1].append(record)
    return tuple(tuple(group) for group in groups)


def _candidate_boundaries(
    records: tuple[_IncidenceDistance, ...],
    attribute: str,
    policy: ShellResolutionPolicy,
) -> tuple[_BoundaryCandidate, ...]:
    groups = _distance_groups(records, attribute, policy.distance_group_tolerance)
    candidates: list[_BoundaryCandidate] = []
    inside: list[_IncidenceDistance] = []
    for index, group in enumerate(groups[:-1]):
        inside.extend(group)
        last_inside = max(getattr(item, attribute) for item in group)
        first_outside = min(getattr(item, attribute) for item in groups[index + 1])
        if last_inside > policy.candidate_rho_max:
            continue
        relative_gap = (first_outside - last_inside) / last_inside
        total_weight = sum(item.incidence.incidence_multiplicity_per_center for item in inside)
        weighted_mean = math.fsum(
            getattr(item, attribute) * item.incidence.incidence_multiplicity_per_center
            for item in inside
        ) / total_weight
        weighted_variance = math.fsum(
            item.incidence.incidence_multiplicity_per_center
            * (getattr(item, attribute) - weighted_mean) ** 2
            for item in inside
        ) / total_weight
        candidates.append(
            _BoundaryCandidate(
                tuple(sorted(item.incidence.incidence_orbit_id for item in inside)),
                relative_gap,
                math.sqrt(weighted_variance) / weighted_mean,
            )
        )
    return tuple(candidates)


def _resolve_projection(
    records: tuple[_IncidenceDistance, ...],
    attribute: str,
    policy: ShellResolutionPolicy,
) -> _BoundaryDecision:
    groups = _distance_groups(records, attribute, policy.distance_group_tolerance)
    if len(groups) < 2:
        diagnostic = Diagnostic(
            Severity.WARNING,
            "crystal_chemistry.shell.candidates_insufficient",
            "At least two observed distance groups are required",
        )
        return _BoundaryDecision(ResolutionStatus.INCOMPLETE, None, (), (diagnostic,))
    candidates = _candidate_boundaries(records, attribute, policy)
    significant = tuple(
        item for item in candidates
        if item.relative_gap >= policy.minimum_shell_gap
    )
    if not significant:
        diagnostic = Diagnostic(
            Severity.WARNING,
            "crystal_chemistry.shell.search_boundary_not_observed",
            "No significant shell boundary was followed by an observed outer group",
        )
        return _BoundaryDecision(ResolutionStatus.INCOMPLETE, None, candidates, (diagnostic,))
    best_gap = max(item.relative_gap for item in significant)
    survivors = tuple(
        item for item in significant
        if best_gap - item.relative_gap <= policy.ambiguity_tolerance
    )
    best_spread = min(item.internal_spread for item in survivors)
    survivors = tuple(
        item for item in survivors
        if item.internal_spread - best_spread <= policy.ambiguity_tolerance
    )
    if len(survivors) == 1:
        return _BoundaryDecision(ResolutionStatus.RESOLVED, survivors[0], significant, ())
    diagnostic = Diagnostic(
        Severity.WARNING,
        "crystal_chemistry.shell.boundary_ambiguous",
        "Distance evidence does not distinguish candidate shell boundaries",
    )
    return _BoundaryDecision(ResolutionStatus.AMBIGUOUS, None, survivors, (diagnostic,))


def _alternative(
    candidate: _BoundaryCandidate,
    records: tuple[_IncidenceDistance, ...],
) -> CoordinationShellAlternative:
    by_id = {item.incidence.incidence_orbit_id: item.incidence for item in records}
    primary = candidate.primary_ids
    secondary = tuple(sorted(set(by_id) - set(primary)))
    geometric_cn = sum(by_id[item].incidence_multiplicity_per_center for item in primary)
    occupied = math.fsum(
        by_id[item].incidence_multiplicity_per_center
        * by_id[item].effective_neighbor_occupancy
        for item in primary
    )
    evidence = (
        SecondaryEvidence(
            "distance_gap",
            EvidenceStatus.SUPPORTIVE,
            f"relative gap {candidate.relative_gap:.12g}; weighted internal spread {candidate.internal_spread:.12g}",
        ),
    )
    alternative_id = "coordination-shell-alternative:" + _digest(
        {"primary_incidence_ids": primary, "secondary_incidence_ids": secondary}
    )
    return CoordinationShellAlternative(
        alternative_id,
        primary,
        secondary,
        geometric_cn,
        occupied,
        evidence,
        ResolutionStatus.RESOLVED,
    )


def _unique_candidates(*groups: tuple[_BoundaryCandidate, ...]) -> tuple[_BoundaryCandidate, ...]:
    by_primary: dict[tuple[str, ...], _BoundaryCandidate] = {}
    for group in groups:
        for candidate in group:
            by_primary.setdefault(candidate.primary_ids, candidate)
    return tuple(by_primary[key] for key in sorted(by_primary))


class CoordinationShellOrbitResolver:
    """Resolve coordination-shell boundaries without materializing contacts."""

    def __init__(self, policy: ShellResolutionPolicy) -> None:
        if not isinstance(policy, ShellResolutionPolicy):
            raise TypeError("policy must be ShellResolutionPolicy")
        self.policy = policy

    def resolve(
        self,
        pair_table: SymmetryPairTable,
        contact_orbits: tuple[ResolvedContactOrbit, ...],
        incidences: tuple[ContactIncidenceOrbit, ...],
    ) -> tuple[CoordinationShellOrbit, ...]:
        if not isinstance(pair_table, SymmetryPairTable):
            raise TypeError("pair_table must be SymmetryPairTable")
        interpretations = {
            item.interpretation_id: item
            for orbit in contact_orbits
            for item in orbit.interpretations
        }
        grouped: dict[tuple[str, str], list[_IncidenceDistance]] = {}
        context_payloads: dict[str, tuple[object, ...]] = {}
        for incidence in incidences:
            try:
                interpretation = interpretations[incidence.interpretation_id]
            except KeyError as exc:
                raise ValueError("incidence references an unknown contact interpretation") from exc
            if interpretation.orientation_mode is not OrientationMode.ENDPOINT_ROLES:
                continue
            context_payload = (
                interpretation.interaction_type.value,
                interpretation.interaction_layer.value,
                interpretation.grammar_priority.value,
            )
            context_id = _interpretation_context_id(interpretation)
            context_payloads[context_id] = context_payload
            if interpretation.normalized_distance_range is None:
                lower = central = upper = math.nan
            else:
                lower, upper = interpretation.normalized_distance_range
                central = math.fsum(
                    item.normalized_distance
                    for item in interpretation.component_pair_interpretations
                ) / len(interpretation.component_pair_interpretations)
            grouped.setdefault((incidence.center_independent_site_id, context_id), []).append(
                _IncidenceDistance(incidence, lower, central, upper)
            )

        shells: list[CoordinationShellOrbit] = []
        for (center_id, context_id), raw_records in sorted(grouped.items()):
            records = tuple(sorted(raw_records, key=lambda item: item.incidence.incidence_orbit_id))
            incomplete = pair_table.status is PairTableStatus.INCOMPLETE or any(
                item.incidence.status is ResolutionStatus.INCOMPLETE
                or not math.isfinite(item.central)
                for item in records
            )
            diagnostics: tuple[Diagnostic, ...]
            selected_id: str | None = None
            if incomplete:
                status = ResolutionStatus.INCOMPLETE
                alternatives: tuple[CoordinationShellAlternative, ...] = ()
                diagnostics = (
                    Diagnostic(
                        Severity.WARNING,
                        "crystal_chemistry.shell.evidence_incomplete",
                        "A complete coordination-shell boundary cannot be derived from the available incidences",
                    ),
                )
            else:
                central = _resolve_projection(records, "central", self.policy)
                lower = _resolve_projection(records, "lower", self.policy)
                upper = _resolve_projection(records, "upper", self.policy)
                bounds_agree = (
                    lower.status is central.status is upper.status
                    and (
                        central.status is not ResolutionStatus.RESOLVED
                        or lower.selected is not None
                        and central.selected is not None
                        and upper.selected is not None
                        and lower.selected.primary_ids
                        == central.selected.primary_ids
                        == upper.selected.primary_ids
                    )
                )
                if not bounds_agree and (
                    lower.status is ResolutionStatus.RESOLVED
                    or upper.status is ResolutionStatus.RESOLVED
                ):
                    candidates = _unique_candidates(
                        (() if lower.selected is None else (lower.selected,)),
                        (() if central.selected is None else (central.selected,)),
                        (() if upper.selected is None else (upper.selected,)),
                        lower.alternatives,
                        central.alternatives,
                        upper.alternatives,
                    )
                    alternatives = tuple(_alternative(item, records) for item in candidates)
                    if len(alternatives) >= 2:
                        status = ResolutionStatus.AMBIGUOUS
                        diagnostics = (
                            Diagnostic(
                                Severity.WARNING,
                                "crystal_chemistry.shell.mixed_occupancy_disagreement",
                                "Valid component interpretations imply different shell boundaries",
                            ),
                        )
                    else:
                        status = ResolutionStatus.INCOMPLETE
                        diagnostics = (
                            Diagnostic(
                                Severity.WARNING,
                                "crystal_chemistry.shell.evidence_incomplete",
                                "Component-distance bounds do not establish two genuine shell alternatives",
                            ),
                        )
                else:
                    status = central.status
                    candidates = central.alternatives
                    alternatives = tuple(_alternative(item, records) for item in candidates)
                    diagnostics = central.diagnostics
                    if central.selected is not None:
                        selected = _alternative(central.selected, records)
                        alternatives_by_id = {item.alternative_id: item for item in alternatives}
                        alternatives_by_id[selected.alternative_id] = selected
                        alternatives = tuple(
                            alternatives_by_id[key] for key in sorted(alternatives_by_id)
                        )
                        selected_id = selected.alternative_id
            shell_id = "coordination-shell-orbit:" + _digest(
                {"center_site_id": center_id, "interpretation_context_id": context_id}
            )
            shells.append(
                CoordinationShellOrbit(
                    shell_id,
                    center_id,
                    context_id,
                    selected_id,
                    alternatives,
                    status,
                    diagnostics,
                    (
                        ("method", "cristma.coordination_shell_orbits:1"),
                        ("interaction_context", context_payloads[context_id]),
                    ),
                )
            )
        return tuple(sorted(shells, key=lambda item: item.shell_orbit_id))


__all__ = [
    "CoordinationShellAlternative",
    "CoordinationShellOrbit",
    "CoordinationShellOrbitResolver",
    "ShellRole",
]
