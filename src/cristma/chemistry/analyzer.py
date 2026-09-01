"""Convenience composition-analysis tool."""

from __future__ import annotations

from dataclasses import dataclass

from cristma.reference_data import ReferenceData

from .classification import ChemicalClassification, classify_composition
from .composition import Composition
from .grammar import CompositionGrammar, compile_composition_grammar


@dataclass(frozen=True, slots=True)
class ChemistryResult:
    composition: Composition
    classification: ChemicalClassification
    grammar: CompositionGrammar


class ChemistryAnalyzer:
    """Stateless tool configured with one immutable reference bundle."""

    def __init__(self, reference: ReferenceData | None = None) -> None:
        self.reference = reference or ReferenceData.default()

    def analyze(self, composition: Composition) -> ChemistryResult:
        classification = classify_composition(composition, self.reference)
        grammar = compile_composition_grammar(composition, classification, self.reference)
        return ChemistryResult(composition, classification, grammar)


__all__ = ["ChemistryAnalyzer", "ChemistryResult"]
