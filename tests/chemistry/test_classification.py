from __future__ import annotations

import math

import pytest

from cristma.chemistry import (
    ChemicalDomain,
    Composition,
    CompositionKind,
    classify_composition,
)
from cristma.reference_data import ReferenceData


@pytest.mark.parametrize(
    ("formula", "kind", "domain", "family"),
    (
        ({"Fe": 1}, CompositionKind.ELEMENTAL, ChemicalDomain.NONE, "elemental.metallic"),
        ({"Si": 1}, CompositionKind.ELEMENTAL, ChemicalDomain.NONE, "elemental.covalent"),
        ({"Fe": 1, "Al": 1}, CompositionKind.COMPOUND, ChemicalDomain.INORGANIC, "inorganic.intermetallic"),
        ({"Ca": 1, "O": 1}, CompositionKind.COMPOUND, ChemicalDomain.INORGANIC, "inorganic.oxide"),
        ({"Na": 1, "Cl": 1}, CompositionKind.COMPOUND, ChemicalDomain.INORGANIC, "inorganic.halide"),
        ({"Fe": 1, "S": 2}, CompositionKind.COMPOUND, ChemicalDomain.INORGANIC, "inorganic.chalcogenide"),
        ({"C": 1, "H": 4}, CompositionKind.COMPOUND, ChemicalDomain.ORGANIC, "organic.molecular"),
        (
            {"Zn": 1, "C": 2, "H": 4, "N": 2},
            CompositionKind.COMPOUND,
            ChemicalDomain.METAL_ORGANIC,
            "coordination.metal_organic_discrete",
        ),
        (
            {"Fe": 1, "C": 5, "H": 5},
            CompositionKind.COMPOUND,
            ChemicalDomain.METAL_ORGANIC,
            "organic.organometallic_molecular",
        ),
    ),
)
def test_ordinary_composition_has_one_actionable_primary_family(
    formula: dict[str, float],
    kind: CompositionKind,
    domain: ChemicalDomain,
    family: str,
) -> None:
    result = classify_composition(Composition.from_mapping(formula), ReferenceData.default())

    assert result.composition_kind is kind
    assert result.domain is domain
    assert result.primary_family == family
    assert result.alternative_families == ()
    assert result.evidence
    assert result.reference_version == "3.1.0-draft"


def test_classification_rejects_invalid_confidence() -> None:
    from cristma.chemistry.classification import ChemicalClassification

    with pytest.raises(ValueError, match="confidence"):
        ChemicalClassification(
            composition_kind=CompositionKind.ELEMENTAL,
            domain=ChemicalDomain.NONE,
            primary_family="elemental.metallic",
            confidence=math.nan,
            evidence=(),
            reference_version="3.1.0-draft",
        )


def test_unknown_inorganic_composition_is_explicitly_unresolved() -> None:
    result = classify_composition(
        Composition.from_mapping({"He": 1, "Ne": 1}),
        ReferenceData.default(),
    )

    assert result.primary_family is None
    assert result.domain is ChemicalDomain.UNRESOLVED
    assert result.diagnostics[0].code == "chemistry.family_unresolved"


def test_oxide_halide_is_routed_as_mixed_anion_not_plain_oxide() -> None:
    result = classify_composition(
        Composition.from_mapping({"La": 1, "O": 1, "F": 1}),
        ReferenceData.default(),
    )

    assert result.primary_family == "inorganic.mixed_anion"
    assert result.alternative_families == ("inorganic.oxide", "inorganic.halide")
