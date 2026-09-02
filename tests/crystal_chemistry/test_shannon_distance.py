from __future__ import annotations

import pytest

from cristma.crystal_chemistry import EvidenceStatus
from cristma.crystal_chemistry.shannon_distance import ShannonDistanceValidator
from cristma.reference_data import ShannonRadii


def _caf2_radii():
    catalog = ShannonRadii.default()
    return (
        catalog.get_exact("Ca", oxidation_state=2, coordination="VIII"),
        catalog.get_exact("F", oxidation_state=-1, coordination="IV"),
    )


def test_shannon_distance_reports_plausible_contact_without_selecting_it() -> None:
    ca, fluorine = _caf2_radii()

    check = ShannonDistanceValidator(minimum_ratio=0.80).evaluate(
        distance=2.36,
        first=ca,
        second=fluorine,
    )

    assert check.status is EvidenceStatus.SUPPORTIVE
    assert check.radius_sum == pytest.approx(ca.ionic_radius + fluorine.ionic_radius)
    assert check.minimum_distance == pytest.approx(0.80 * check.radius_sum)
    assert check.distance_ratio == pytest.approx(2.36 / check.radius_sum)
    assert check.excludes_contact is False


def test_shannon_distance_marks_only_explicit_overlap_as_contradictory() -> None:
    ca, fluorine = _caf2_radii()

    check = ShannonDistanceValidator(minimum_ratio=0.90).evaluate(
        distance=1.80,
        first=ca,
        second=fluorine,
    )

    assert check.status is EvidenceStatus.CONTRADICTORY
    assert check.excludes_contact is False


def test_shannon_distance_policy_is_explicit_and_validated() -> None:
    validator = ShannonDistanceValidator(minimum_ratio=0.80)

    assert validator.get_config() == {"minimum_ratio": 0.80}
    assert validator.clone(minimum_ratio=0.85).minimum_ratio == pytest.approx(0.85)
    assert validator.minimum_ratio == pytest.approx(0.80)
    with pytest.raises(ValueError, match="minimum_ratio"):
        ShannonDistanceValidator(minimum_ratio=0.0)
