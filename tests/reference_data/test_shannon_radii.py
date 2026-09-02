from __future__ import annotations

import pytest

from cristma.reference_data import ReferenceData, ShannonRadii, ShannonSpinState


def test_shannon_lookup_distinguishes_spin_states_without_guessing() -> None:
    radii = ShannonRadii.default()

    alternatives = radii.find("Fe", oxidation_state=2, coordination="VI")

    assert tuple(item.spin_state for item in alternatives) == (
        ShannonSpinState.HIGH_SPIN,
        ShannonSpinState.LOW_SPIN,
    )
    assert tuple(item.ionic_radius for item in alternatives) == pytest.approx((0.78, 0.61))
    assert radii.get_exact(
        "Fe",
        oxidation_state=2,
        coordination="VI",
        spin_state=ShannonSpinState.HIGH_SPIN,
    ).crystal_radius == pytest.approx(0.92)


def test_shannon_lookup_preserves_coordination_geometry_labels() -> None:
    square_planar = ShannonRadii.default().get_exact(
        "Ag", oxidation_state=1, coordination="IVSQ"
    )

    assert square_planar.coordination == "IVSQ"
    assert square_planar.ionic_radius == pytest.approx(1.02)


def test_shannon_exact_lookup_never_selects_an_ambiguous_record() -> None:
    radii = ShannonRadii.default()

    with pytest.raises(LookupError, match="2 Shannon radii"):
        radii.get_exact("Fe", oxidation_state=2, coordination="VI")
    assert radii.find("He", oxidation_state=2, coordination="VI") == ()


def test_reference_data_exposes_shannon_catalog() -> None:
    radius = ReferenceData.default().shannon_radii.get_exact(
        "Ca", oxidation_state=2, coordination="VIII"
    )

    assert radius.ionic_radius == pytest.approx(1.12)
    assert radius.unit == "angstrom"
    assert radius.dataset_id == "cristma.shannon_radii.pymatgen"


def test_shannon_catalog_has_complete_unique_pinned_dataset() -> None:
    records = ShannonRadii.default().records

    assert len(records) == 493
    assert len({
        (record.symbol, record.oxidation_state, record.coordination, record.spin_state)
        for record in records
    }) == 493
