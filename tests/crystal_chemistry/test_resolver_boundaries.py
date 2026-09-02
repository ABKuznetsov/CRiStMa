from __future__ import annotations

import pytest

from cristma.crystal_chemistry import EvidenceStatus, ResolutionStatus, ShellResolutionPolicy
from cristma.crystal_chemistry.resolver import _resolve_rho_values


POLICY = ShellResolutionPolicy(1.6, 0.01, 0.08, 0.01, 2.0)


def codes(decision) -> set[str]:
    return {item.code for item in decision.diagnostics}


def test_candidate_uses_relative_gap_and_internal_spread() -> None:
    decision = _resolve_rho_values((0.98, 1.00, 1.01, 1.28), POLICY)

    assert decision.status is ResolutionStatus.RESOLVED
    assert decision.selected is not None
    assert decision.selected.geometric_CN == 3
    assert decision.selected.relative_gap == pytest.approx((1.28 - 1.01) / 1.01)
    assert decision.selected.internal_spread == pytest.approx((1.01 - 0.98) / 1.00)


def test_close_non_dominated_boundaries_remain_ambiguous() -> None:
    tolerant = ShellResolutionPolicy(1.6, 0.01, 0.08, 0.25, 2.0)

    decision = _resolve_rho_values((0.95, 1.00, 1.10, 1.16, 1.28), tolerant)

    assert decision.status is ResolutionStatus.AMBIGUOUS
    assert decision.selected is None
    assert len(decision.alternatives) >= 2
    assert "crystal_chemistry.shell.boundary_ambiguous" in codes(decision)


def test_candidate_range_without_significant_outer_gap_is_incomplete() -> None:
    decision = _resolve_rho_values((0.98, 1.00, 1.01), POLICY)

    assert decision.status is ResolutionStatus.INCOMPLETE
    assert codes(decision) == {
        "crystal_chemistry.shell.search_boundary_not_observed"
    }


def test_single_distance_group_is_insufficient() -> None:
    decision = _resolve_rho_values((1.00,), POLICY)

    assert decision.status is ResolutionStatus.INCOMPLETE
    assert codes(decision) == {"crystal_chemistry.shell.candidates_insufficient"}


def test_secondary_evidence_is_explicit_and_does_not_select_boundary() -> None:
    decision = _resolve_rho_values((0.98, 1.00, 1.01, 1.28), POLICY)

    assert {(item.method, item.status) for item in decision.evidence} == {
        ("bvs", EvidenceStatus.NOT_AVAILABLE),
        ("coordination_geometry", EvidenceStatus.NOT_APPLICABLE),
    }


def test_close_values_form_one_distance_group() -> None:
    decision = _resolve_rho_values((1.000, 1.004, 1.009, 1.20), POLICY)

    assert decision.status is ResolutionStatus.RESOLVED
    assert decision.selected is not None
    assert decision.selected.boundary_group == 0
    assert decision.selected.geometric_CN == 3
