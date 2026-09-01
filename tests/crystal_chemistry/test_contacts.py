from __future__ import annotations

import math

import pytest

from cristma.chemistry import GrammarOperation, InteractionPriority
from cristma.core.values import MeasuredValue
from cristma.crystallography import GeometricContact
from cristma.crystal_chemistry import (
    ComponentPairInterpretation,
    ContactClassification,
    CoordinationShell,
    EvidenceStatus,
    ResolutionStatus,
    ResolvedContact,
    SecondaryEvidence,
    ShellResolutionPolicy,
)
from cristma.structure import SiteComponent


def interpretation() -> ComponentPairInterpretation:
    ca = SiteComponent("Ca", MeasuredValue(1.0, None, "1"))
    oxygen = SiteComponent("O", MeasuredValue(0.75, None, "0.75"))
    return ComponentPairInterpretation(
        first_species=ca.species,
        second_species=oxygen.species,
        first_occupancy=1.0,
        second_occupancy=0.75,
        radius_sum=2.42,
        normalized_distance=1.0,
        occupancy_weight=0.75,
        interaction_type=GrammarOperation.CENTRE_LIGAND_SHELL,
        grammar_priority=InteractionPriority.PRIMARY,
    )


def resolved_contact(index: int) -> ResolvedContact:
    pair = interpretation()
    geometric = GeometricContact(
        contact_id=f"contact:M|O{index}|finite",
        first_atom_id="M",
        second_atom_id=f"O{index}",
        cell_translation=None,
        distance=2.42,
        vector_cartesian=(2.42, 0.0, 0.0),
        first_source_site_id="site:M",
        second_source_site_id=f"site:O{index}",
        geometric_provenance="neighbor_graph",
    )
    return ResolvedContact(
        geometric_contact=geometric,
        interaction_type=GrammarOperation.CENTRE_LIGAND_SHELL,
        grammar_priority=InteractionPriority.PRIMARY,
        contact_classification=ContactClassification.PRIMARY,
        component_interpretations=(pair,),
        normalized_distance_min=1.0,
        normalized_distance_max=1.0,
        neighbor_total_occupancy=0.75,
        evidence=(SecondaryEvidence("bvs", EvidenceStatus.NOT_AVAILABLE, "not run"),),
        provenance=(("method", "fixture"),),
    )


def test_policy_is_explicit_cloneable_and_dimensionless() -> None:
    policy = ShellResolutionPolicy(1.45, 0.01, 0.08, 0.01)

    clone = policy.clone(candidate_rho_max=1.60)

    assert policy.get_config() == {
        "candidate_rho_max": 1.45,
        "distance_group_tolerance": 0.01,
        "minimum_shell_gap": 0.08,
        "ambiguity_tolerance": 0.01,
    }
    assert clone.candidate_rho_max == 1.60
    assert policy.candidate_rho_max == 1.45


@pytest.mark.parametrize("value", [0.0, -0.1, math.inf, math.nan, True])
def test_policy_rejects_nonpositive_nonfinite_and_boolean_values(value) -> None:
    with pytest.raises(ValueError):
        ShellResolutionPolicy(value, 0.01, 0.08, 0.01)


def test_shell_counts_geometric_positions_and_occupancy_separately() -> None:
    contacts = tuple(resolved_contact(index) for index in range(4))

    shell = CoordinationShell.resolved("site:M", "M", contacts)

    assert shell.status is ResolutionStatus.RESOLVED
    assert shell.geometric_CN == 4
    assert shell.mean_occupied_neighbors == pytest.approx(3.0)


def test_resolved_contact_bounds_must_cover_every_interpretation() -> None:
    contact = resolved_contact(1)

    with pytest.raises(ValueError, match="normalized-distance bounds"):
        ResolvedContact(
            geometric_contact=contact.geometric_contact,
            interaction_type=contact.interaction_type,
            grammar_priority=contact.grammar_priority,
            contact_classification=contact.contact_classification,
            component_interpretations=contact.component_interpretations,
            normalized_distance_min=0.8,
            normalized_distance_max=0.9,
            neighbor_total_occupancy=0.75,
            evidence=contact.evidence,
            provenance=contact.provenance,
        )


def test_component_interpretation_keeps_interaction_context() -> None:
    pair = interpretation()

    assert pair.species_symbols == ("Ca", "O")
    assert pair.interaction_type is GrammarOperation.CENTRE_LIGAND_SHELL
    assert pair.grammar_priority is InteractionPriority.PRIMARY
