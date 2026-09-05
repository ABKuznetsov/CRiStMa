from __future__ import annotations

from dataclasses import replace
from types import MappingProxyType

import pytest

from cristma.reference_data import (
    ReferenceData,
    load_chemical_reference,
    validate_reference_integrity,
)


def test_default_reference_is_validated_v31() -> None:
    reference = ReferenceData.default()

    assert reference.chemical.schema_name == "CrIStMa Chemical Reference DB"
    assert reference.chemical.schema_version == "3.1.0-draft"
    assert reference.chemical.family("inorganic.oxide")["profile_id"] == "inorganic.oxide"
    assert reference.chemical.boundary_case("CaSi2")["refined"]["preferred_candidates"] == (
        "inorganic.tetrelide",
        "inorganic.zintl",
    )


def test_historical_v3_counts_are_stable() -> None:
    report = validate_reference_integrity(load_chemical_reference("3.0.0-draft"))

    assert (report.family_count, report.group_count, report.boundary_case_count) == (
        103,
        243,
        155,
    )
    assert report.valid


def test_reference_data_is_deeply_immutable() -> None:
    family = ReferenceData.default().chemical.family("inorganic.carbide.covalent")

    with pytest.raises(TypeError):
        family["profile_id"] = "changed"


def test_reference_rejects_unknown_version() -> None:
    with pytest.raises(ValueError, match="unsupported chemical reference version"):
        load_chemical_reference("4")


def test_v31_routes_concrete_grammar_templates() -> None:
    reference = ReferenceData.default().chemical

    assert reference.grammar_route("inorganic.halide.bromide") == "ionic_halide"
    assert reference.composition_family_route("halide") == "inorganic.halide"
    assert reference.grammar_template("ionic_halide")["interactions"][0]["operation"] == "centre_ligand_shell"


def test_reference_rejects_unknown_interaction_layer() -> None:
    reference = ReferenceData.default().chemical
    templates = dict(reference._grammar_templates)
    template = dict(templates["ionic_halide"])
    interactions = list(template["interactions"])
    interactions[0] = MappingProxyType({**interactions[0], "layer": "unknown"})
    template["interactions"] = tuple(interactions)
    templates["ionic_halide"] = MappingProxyType(template)
    invalid = replace(reference, _grammar_templates=MappingProxyType(templates))

    with pytest.raises(ValueError, match="unknown grammar layer"):
        validate_reference_integrity(invalid)
