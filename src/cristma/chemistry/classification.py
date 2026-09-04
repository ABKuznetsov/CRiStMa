"""Deterministic composition-level material classification."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
import math

from cristma.diagnostics import Diagnostic, Severity
from cristma.reference_data import ReferenceData

from .composition import Composition
from .evidence import ChemicalEvidence


class CompositionKind(StrEnum):
    ELEMENTAL = "elemental"
    COMPOUND = "compound"


class ChemicalDomain(StrEnum):
    NONE = "none"
    INORGANIC = "inorganic"
    ORGANIC = "organic"
    METAL_ORGANIC = "metal_organic"
    UNRESOLVED = "unresolved"


@dataclass(frozen=True, slots=True)
class ChemicalClassification:
    composition_kind: CompositionKind
    domain: ChemicalDomain
    primary_family: str | None
    confidence: float
    evidence: tuple[ChemicalEvidence, ...]
    reference_version: str
    alternative_families: tuple[str, ...] = ()
    modifiers: tuple[str, ...] = ()
    diagnostics: tuple[Diagnostic, ...] = ()
    method_id: str = "cristma.composition_classification"
    method_version: str = "1"

    def __post_init__(self) -> None:
        if not math.isfinite(self.confidence) or not 0 <= self.confidence <= 1:
            raise ValueError("classification confidence must be finite and between zero and one")
        if not self.evidence:
            raise ValueError("classification requires evidence")
        if self.composition_kind is CompositionKind.ELEMENTAL and self.domain is not ChemicalDomain.NONE:
            raise ValueError("elemental classification requires domain NONE")
        if self.composition_kind is CompositionKind.COMPOUND and self.domain is ChemicalDomain.NONE:
            raise ValueError("compound classification cannot use domain NONE")
        if self.primary_family is not None and self.primary_family in self.alternative_families:
            raise ValueError("primary family cannot also be an alternative")


_ELEMENTAL_COVALENT = frozenset({"B", "C", "Si", "Ge"})
_ELEMENTAL_MOLECULAR = frozenset({"H", "N", "O", "F", "P", "S", "Cl", "Se", "Br", "I"})
_ORGANIC_DONORS = frozenset({"N", "O", "S", "P", "F", "Cl", "Br", "I"})


def _evidence(code: str, message: str, elements: frozenset[str]) -> tuple[ChemicalEvidence, ...]:
    return (ChemicalEvidence(code, message, tuple(sorted(elements))),)


def _result(
    *,
    elements: frozenset[str],
    kind: CompositionKind,
    domain: ChemicalDomain,
    family: str | None,
    reference: ReferenceData,
    code: str,
    message: str,
    confidence: float = 0.85,
    alternatives: tuple[str, ...] = (),
    diagnostics: tuple[Diagnostic, ...] = (),
) -> ChemicalClassification:
    if family is not None:
        reference.chemical.family(family)
    for alternative in alternatives:
        reference.chemical.family(alternative)
    return ChemicalClassification(
        composition_kind=kind,
        domain=domain,
        primary_family=family,
        alternative_families=alternatives,
        confidence=confidence,
        evidence=_evidence(code, message, elements),
        diagnostics=diagnostics,
        reference_version=reference.chemical.schema_version,
    )


def classify_composition(
    composition: Composition,
    reference: ReferenceData | None = None,
) -> ChemicalClassification:
    """Return one actionable family for ordinary compositions."""

    reference = reference or ReferenceData.default()
    elements = frozenset(composition.elements)
    metals = frozenset(symbol for symbol in elements if reference.elements.by_symbol(symbol).is_metal)

    if len(elements) == 1:
        symbol = next(iter(elements))
        if reference.elements.by_symbol(symbol).is_metal:
            family = "elemental.metallic"
        elif symbol in _ELEMENTAL_MOLECULAR:
            family = "elemental.molecular"
        else:
            family = "elemental.covalent"
        return _result(
            elements=elements,
            kind=CompositionKind.ELEMENTAL,
            domain=ChemicalDomain.NONE,
            family=family,
            reference=reference,
            code="chemistry.elemental_family",
            message=f"Elemental {symbol} is assigned to {family}.",
        )

    if {"C", "H"} <= elements:
        if metals:
            family = (
                "coordination.metal_organic_discrete"
                if elements & _ORGANIC_DONORS
                else "organic.organometallic_molecular"
            )
            return _result(
                elements=elements,
                kind=CompositionKind.COMPOUND,
                domain=ChemicalDomain.METAL_ORGANIC,
                family=family,
                reference=reference,
                code="chemistry.metal_organic_composition",
                message="Metal plus carbon-hydrogen composition selects metal-organic analysis.",
                confidence=0.7,
            )
        return _result(
            elements=elements,
            kind=CompositionKind.COMPOUND,
            domain=ChemicalDomain.ORGANIC,
            family="organic.molecular",
            reference=reference,
            code="chemistry.organic_composition",
            message="Carbon-hydrogen composition selects molecular organic analysis.",
            confidence=0.8,
        )

    if metals == elements:
        return _result(
            elements=elements,
            kind=CompositionKind.COMPOUND,
            domain=ChemicalDomain.INORGANIC,
            family="inorganic.intermetallic",
            reference=reference,
            code="chemistry.intermetallic_composition",
            message="All occupied elements are metallic; use intermetallic coordination.",
        )

    halogens = reference.chemical.element_set("halogens")
    chalcogens = reference.chemical.element_set("chalcogens") - {"O"}
    pnictogens = reference.chemical.element_set("pnictogens") - {"N"}
    family: str | None = None
    code = "chemistry.inorganic_family"
    secondary_anion_families: list[str] = []
    if elements & halogens:
        secondary_anion_families.append("inorganic.halide")
    if elements & chalcogens:
        secondary_anion_families.append("inorganic.chalcogenide")
    if "N" in elements:
        secondary_anion_families.append("inorganic.nitride")
    if elements & pnictogens:
        secondary_anion_families.append("inorganic.pnictide")

    if "O" in elements and secondary_anion_families:
        return _result(
            elements=elements,
            kind=CompositionKind.COMPOUND,
            domain=ChemicalDomain.INORGANIC,
            family="inorganic.mixed_anion",
            alternatives=("inorganic.oxide", *tuple(dict.fromkeys(secondary_anion_families))),
            reference=reference,
            code="chemistry.mixed_anion_candidate",
            message="Oxygen and another anion-forming subsystem select mixed-anion analysis.",
            confidence=0.75,
        )

    if "O" in elements:
        family = "inorganic.oxide"
    elif elements & halogens:
        family = "inorganic.halide"
    elif "N" in elements:
        family = "inorganic.nitride"
    elif elements & chalcogens:
        family = "inorganic.chalcogenide"
    elif elements & pnictogens:
        family = "inorganic.pnictide"
    elif "C" in elements:
        family = "inorganic.carbide"
    elif "B" in elements:
        family = "inorganic.boride"
    elif elements & {"Si", "Ge", "Sn", "Pb"}:
        family = "inorganic.tetrelide"

    if family is not None:
        return _result(
            elements=elements,
            kind=CompositionKind.COMPOUND,
            domain=ChemicalDomain.INORGANIC,
            family=family,
            reference=reference,
            code=code,
            message=f"Composition selects {family} analysis.",
        )

    return _result(
        elements=elements,
        kind=CompositionKind.COMPOUND,
        domain=ChemicalDomain.UNRESOLVED,
        family=None,
        reference=reference,
        code="chemistry.family_unresolved",
        message="The implemented composition rules do not identify this material family.",
        confidence=0.0,
        diagnostics=(
            Diagnostic(
                Severity.WARNING,
                "chemistry.family_unresolved",
                "No composition-level material family was identified.",
            ),
        ),
    )


__all__ = [
    "ChemicalClassification",
    "ChemicalDomain",
    "CompositionKind",
    "classify_composition",
]
