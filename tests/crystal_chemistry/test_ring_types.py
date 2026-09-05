from __future__ import annotations

import pytest

from cristma.chemistry import Composition
from cristma.crystal_chemistry.rings import (
    PeriodicUnitRef,
    RingAnalysisResult,
    RingAnalysisStatus,
    RingSearchPolicy,
    StructuralRing,
    StructuralRingOrbit,
)
from cristma.diagnostics import Diagnostic, Severity
from cristma.structure import PeriodicAtomRef


def _ring(ring_id: str = "ring:abc") -> StructuralRing:
    return StructuralRing(
        ring_id=ring_id,
        parent_block_id="block:framework",
        representation_id="representation:structural",
        unit_refs=(
            PeriodicUnitRef("unit:A", (0, 0, 0)),
            PeriodicUnitRef("unit:B", (0, 0, 0)),
            PeriodicUnitRef("unit:C", (0, 0, 0)),
        ),
        connection_ids=("edge:AB", "edge:BC", "edge:CA"),
        connector_atom_refs=(PeriodicAtomRef("atom:O1", (0, 0, 0)),),
        composition=Composition.from_mapping({"B": 3, "O": 7}),
        translation_sum=(0, 0, 0),
        provenance=(("method", "test"),),
    )


def test_structural_ring_exposes_its_cycle_size() -> None:
    assert _ring().size == 3


def test_structural_ring_rejects_a_winding_cycle() -> None:
    with pytest.raises(ValueError, match="zero translation"):
        StructuralRing(
            **{
                **{field: getattr(_ring(), field) for field in (
                    "ring_id",
                    "parent_block_id",
                    "representation_id",
                    "unit_refs",
                    "connection_ids",
                    "connector_atom_refs",
                    "composition",
                    "provenance",
                )},
                "translation_sum": (1, 0, 0),
            }
        )


def test_ring_policy_is_cloneable_and_rejects_invalid_limits() -> None:
    policy = RingSearchPolicy(maximum_ring_size=12, maximum_states_per_connection=50_000)

    assert policy.clone(maximum_ring_size=8).maximum_ring_size == 8
    assert policy.get_config()["maximum_ring_size"] == 12
    with pytest.raises(ValueError, match="positive"):
        RingSearchPolicy(maximum_paths_per_connection=0)


def test_incomplete_result_requires_a_diagnostic() -> None:
    with pytest.raises(ValueError, match="diagnostic"):
        RingAnalysisResult((), (), RingAnalysisStatus.INCOMPLETE)


def test_orbit_members_must_exist_in_result_and_share_context() -> None:
    ring = _ring()
    missing = StructuralRingOrbit(
        orbit_id="ring-orbit:missing",
        parent_block_id=ring.parent_block_id,
        representation_id=ring.representation_id,
        representative_ring_id="ring:absent",
        ring_ids=("ring:absent",),
        multiplicity=1,
        composition=ring.composition,
        size=ring.size,
    )

    with pytest.raises(ValueError, match="unknown ring"):
        RingAnalysisResult((ring,), (missing,), RingAnalysisStatus.COMPLETE)

    diagnostic = Diagnostic(
        Severity.WARNING,
        "crystal_chemistry.rings.search_limit_reached",
        "Search limit reached.",
    )
    result = RingAnalysisResult(
        (ring,),
        (),
        RingAnalysisStatus.INCOMPLETE,
        (diagnostic,),
    )
    assert result.diagnostics == (diagnostic,)
