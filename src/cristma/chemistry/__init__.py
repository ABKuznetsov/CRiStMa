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

__all__ = [
    "ChargedSpecies",
    "ChemicalClassification",
    "ChemicalDomain",
    "ChemicalEvidence",
    "ChemicalSpecies",
    "Composition",
    "CompositionKind",
    "ELEMENT_SYMBOLS",
    "ELEMENT_SYMBOLS_BY_ATOMIC_NUMBER",
    "ElementSpecies",
    "IsotopeSpecies",
    "UnknownSpecies",
    "as_species",
    "classify_composition",
    "element_from_atomic_number",
    "normalize_element",
]
