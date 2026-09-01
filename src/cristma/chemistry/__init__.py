"""Chemical identities and crystal-chemistry data."""

from .elements import (
    ELEMENT_SYMBOLS,
    ELEMENT_SYMBOLS_BY_ATOMIC_NUMBER,
    element_from_atomic_number,
    normalize_element,
)
from .species import (
    ChargedSpecies,
    ChemicalSpecies,
    ElementSpecies,
    IsotopeSpecies,
    UnknownSpecies,
    as_species,
)
from .composition import Composition
from .classification import (
    ChemicalClassification,
    ChemicalDomain,
    CompositionKind,
    classify_composition,
)
from .evidence import ChemicalEvidence
from .grammar import (
    CandidateInteraction,
    CompositionGrammar,
    DecompositionMode,
    GrammarOperation,
    InteractionLayer,
    InteractionPriority,
    compile_composition_grammar,
)
from .analyzer import ChemistryAnalyzer, ChemistryResult

__all__ = [
    "CandidateInteraction",
    "ChargedSpecies",
    "ChemicalClassification",
    "ChemicalDomain",
    "ChemicalEvidence",
    "ChemicalSpecies",
    "ChemistryAnalyzer",
    "ChemistryResult",
    "Composition",
    "CompositionKind",
    "CompositionGrammar",
    "DecompositionMode",
    "ELEMENT_SYMBOLS",
    "ELEMENT_SYMBOLS_BY_ATOMIC_NUMBER",
    "ElementSpecies",
    "GrammarOperation",
    "InteractionLayer",
    "InteractionPriority",
    "IsotopeSpecies",
    "UnknownSpecies",
    "as_species",
    "classify_composition",
    "compile_composition_grammar",
    "element_from_atomic_number",
    "normalize_element",
]
