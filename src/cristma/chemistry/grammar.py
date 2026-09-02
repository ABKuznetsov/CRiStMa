"""Concrete composition-level interaction requests for geometry tools."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
import math

from cristma.diagnostics import Diagnostic, Severity
from cristma.reference_data import ReferenceData

from .classification import ChemicalClassification
from .composition import Composition
from .evidence import ChemicalEvidence


class DecompositionMode(StrEnum):
    IONIC_SUBLATTICES = "ionic_sublattices"
    STRUCTURAL_ANION_SUBSYSTEM = "structural_anion_subsystem"
    CATION_ANION_SUBSYSTEM = "cation_anion_subsystem"
    COVALENT_NETWORK = "covalent_network"
    MOLECULAR_COMPONENTS = "molecular_components"
    METAL_ORGANIC = "metal_organic"
    METALLIC_SUBLATTICES = "metallic_sublattices"
    UNRESOLVED = "unresolved"


class GrammarOperation(StrEnum):
    CENTRE_LIGAND_SHELL = "centre_ligand_shell"
    INTERSTITIAL_COORDINATION = "interstitial_coordination"
    MIXED_ANION_COORDINATION = "mixed_anion_coordination"
    COVALENT_NETWORK = "covalent_network"
    INTRA_SUBSYSTEM_BONDS = "intra_subsystem_bonds"
    METALLIC_COORDINATION = "metallic_coordination"


class InteractionLayer(StrEnum):
    STRUCTURAL = "structural"
    INTERSTITIAL = "interstitial"
    COORDINATION = "coordination"
    INTRA_SUBSYSTEM = "intra_subsystem"
    INTRAMOLECULAR = "intramolecular"
    METALLIC = "metallic"


class InteractionPriority(StrEnum):
    PRIMARY = "primary"
    ALLOWED = "allowed"


def _symbols(values: tuple[str, ...]) -> tuple[str, ...]:
    result = tuple(sorted(set(values)))
    if not result:
        raise ValueError("interaction element group must not be empty")
    return result


@dataclass(frozen=True, slots=True)
class CandidateInteraction:
    first_elements: tuple[str, ...]
    second_elements: tuple[str, ...]
    operation: GrammarOperation
    layer: InteractionLayer
    priority: InteractionPriority
    centre_elements: tuple[str, ...]
    ligand_elements: tuple[str, ...]
    evidence: tuple[ChemicalEvidence, ...]

    def __post_init__(self) -> None:
        first = _symbols(self.first_elements)
        second = _symbols(self.second_elements)
        first, second = min((first, second), (second, first))
        centres = _symbols(self.centre_elements)
        ligands = _symbols(self.ligand_elements)
        if not set(centres) <= set(first) | set(second):
            raise ValueError("centre elements must belong to the interaction pair")
        if not set(ligands) <= set(first) | set(second):
            raise ValueError("ligand elements must belong to the interaction pair")
        if not self.evidence:
            raise ValueError("candidate interaction requires evidence")
        object.__setattr__(self, "first_elements", first)
        object.__setattr__(self, "second_elements", second)
        object.__setattr__(self, "centre_elements", centres)
        object.__setattr__(self, "ligand_elements", ligands)


@dataclass(frozen=True, slots=True)
class CompositionGrammar:
    mode: DecompositionMode
    candidate_interactions: tuple[CandidateInteraction, ...]
    confidence: float
    evidence: tuple[ChemicalEvidence, ...]
    diagnostics: tuple[Diagnostic, ...]
    reference_version: str
    method_id: str = "cristma.composition_grammar"
    method_version: str = "1"

    def __post_init__(self) -> None:
        if not math.isfinite(self.confidence) or not 0 <= self.confidence <= 1:
            raise ValueError("grammar confidence must be finite and between zero and one")
        if not self.evidence:
            raise ValueError("composition grammar requires evidence")
        if len(set(self.candidate_interactions)) != len(self.candidate_interactions):
            raise ValueError("composition grammar contains duplicate interactions")


def _interaction(
    first: tuple[str, ...],
    second: tuple[str, ...],
    operation: GrammarOperation,
    priority: InteractionPriority,
    evidence: tuple[ChemicalEvidence, ...],
    *,
    centres: tuple[str, ...] | None = None,
    ligands: tuple[str, ...] | None = None,
    layer: InteractionLayer = InteractionLayer.COORDINATION,
) -> CandidateInteraction:
    return CandidateInteraction(
        first,
        second,
        operation,
        layer,
        priority,
        centres or first,
        ligands or second,
        evidence,
    )


def compile_composition_grammar(
    composition: Composition,
    classification: ChemicalClassification,
    reference: ReferenceData | None = None,
) -> CompositionGrammar:
    """Compile declarative reference-data templates into concrete searches."""

    reference = reference or ReferenceData.default()
    elements = frozenset(composition.elements)
    family = classification.primary_family
    evidence = (
        ChemicalEvidence(
            "chemistry.interaction_grammar",
            f"Interaction search compiled for {family or 'unresolved composition'}.",
            tuple(sorted(elements)),
            (family,) if family else (),
        ),
    )
    if family is None:
        return CompositionGrammar(
            mode=DecompositionMode.UNRESOLVED,
            candidate_interactions=(),
            confidence=classification.confidence,
            evidence=evidence,
            diagnostics=classification.diagnostics,
            reference_version=reference.chemical.schema_version,
            method_version="2",
        )

    route_key = family.rsplit(".", 1)[-1]
    try:
        routed_family = reference.chemical.composition_family_route(
            route_key, elements
        )
        template_id = reference.chemical.grammar_route(routed_family)
    except KeyError:
        try:
            template_id = reference.chemical.grammar_route(family)
        except KeyError:
            diagnostic = Diagnostic(
                Severity.WARNING,
                "chemistry.grammar_unresolved",
                f"No declarative composition grammar is routed for {family}.",
            )
            return CompositionGrammar(
                mode=DecompositionMode.UNRESOLVED,
                candidate_interactions=(),
                confidence=classification.confidence,
                evidence=evidence,
                diagnostics=classification.diagnostics + (diagnostic,),
                reference_version=reference.chemical.schema_version,
                method_version="2",
            )

    template = reference.chemical.grammar_template(template_id)
    records = tuple(template["interactions"])
    selector_names = {
        str(record[field])
        for record in records
        for field in ("first", "second", "centres", "ligands")
    }
    fixed: dict[str, frozenset[str]] = {}
    for selector in selector_names:
        if selector in {
            "all_elements",
            "remaining_elements",
            "remaining_electropositive_elements",
            "metal_elements",
            "nonmetal_elements",
        }:
            continue
        try:
            selected = reference.chemical.grammar_element_set(selector)
        except KeyError:
            selected = reference.chemical.element_set(selector)
        fixed[selector] = elements & selected
    explicitly_selected = frozenset().union(*fixed.values()) if fixed else frozenset()

    def select(identifier: str) -> frozenset[str]:
        if identifier == "all_elements":
            return elements
        if identifier == "remaining_elements":
            return elements - explicitly_selected
        if identifier == "remaining_electropositive_elements":
            return frozenset(
                symbol
                for symbol in elements - explicitly_selected
                if reference.elements.by_symbol(symbol).is_metal or symbol == "H"
            )
        if identifier == "metal_elements":
            return frozenset(
                symbol for symbol in elements
                if reference.elements.by_symbol(symbol).is_metal
            )
        if identifier == "nonmetal_elements":
            return frozenset(
                symbol for symbol in elements
                if not reference.elements.by_symbol(symbol).is_metal
            )
        return fixed[identifier]

    interactions: list[CandidateInteraction] = []
    for record in records:
        first = tuple(sorted(select(str(record["first"]))))
        second = tuple(sorted(select(str(record["second"]))))
        centres = tuple(sorted(select(str(record["centres"]))))
        ligands = tuple(sorted(select(str(record["ligands"]))))
        if not first or not second or not centres or not ligands:
            continue
        operation = GrammarOperation(str(record["operation"]))
        layer = InteractionLayer(str(record["layer"]))
        priority = InteractionPriority(str(record["priority"]))
        if record.get("expand") == "distinct_unordered_pairs":
            ordered = tuple(sorted(set(first) | set(second)))
            for index, first_symbol in enumerate(ordered):
                for second_symbol in ordered[index + 1:]:
                    interactions.append(_interaction(
                        (first_symbol,), (second_symbol,), operation, priority,
                        evidence, centres=(first_symbol,), ligands=(second_symbol,),
                        layer=layer,
                    ))
            continue
        interactions.append(_interaction(
            first, second, operation, priority, evidence,
            centres=centres, ligands=ligands, layer=layer,
        ))

    return CompositionGrammar(
        mode=DecompositionMode(str(template["mode"])),
        candidate_interactions=tuple(interactions),
        confidence=classification.confidence,
        evidence=evidence,
        diagnostics=classification.diagnostics,
        reference_version=reference.chemical.schema_version,
        method_version="2",
    )


__all__ = [
    "CandidateInteraction",
    "CompositionGrammar",
    "DecompositionMode",
    "GrammarOperation",
    "InteractionLayer",
    "InteractionPriority",
    "compile_composition_grammar",
]
