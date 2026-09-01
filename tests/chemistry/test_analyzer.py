from __future__ import annotations

from cristma.chemistry import ChemistryAnalyzer, Composition
from cristma.reference_data import ReferenceData


def test_analyzer_returns_composition_classification_and_grammar() -> None:
    composition = Composition.from_mapping({"Ca": 1, "O": 1})
    analyzer = ChemistryAnalyzer()

    result = analyzer.analyze(composition)

    assert result.composition is composition
    assert result.classification.primary_family == "inorganic.oxide"
    assert result.grammar.candidate_interactions


def test_analyzer_stores_reference_configuration_not_last_result() -> None:
    analyzer = ChemistryAnalyzer(reference=ReferenceData.default())

    analyzer.analyze(Composition.from_mapping({"Fe": 1}))

    assert not hasattr(analyzer, "last_result")
    assert not hasattr(analyzer, "current_structure")
