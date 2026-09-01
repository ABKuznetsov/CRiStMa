"""Concrete composition-level interaction requests for geometry tools."""

from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum
import math

from cristma.diagnostics import Diagnostic
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
    PRIMARY_COORDINATION = "primary_coordination"
    INTRA_SUBSYSTEM = "intra_subsystem"
    INTRAMOLECULAR = "intramolecular"


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
    layer: InteractionLayer = InteractionLayer.PRIMARY_COORDINATION,
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
    """Translate one material family into concrete element-pair searches."""

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
    interactions: list[CandidateInteraction] = []
    mode = DecompositionMode.UNRESOLVED

    if family == "elemental.metallic":
        group = tuple(elements)
        mode = DecompositionMode.METALLIC_SUBLATTICES
        interactions.append(_interaction(group, group, GrammarOperation.METALLIC_COORDINATION, InteractionPriority.PRIMARY, evidence))
    elif family in {"elemental.covalent", "elemental.molecular"}:
        group = tuple(elements)
        mode = DecompositionMode.COVALENT_NETWORK if family.endswith("covalent") else DecompositionMode.MOLECULAR_COMPONENTS
        interactions.append(_interaction(group, group, GrammarOperation.COVALENT_NETWORK, InteractionPriority.PRIMARY, evidence, layer=InteractionLayer.INTRAMOLECULAR))
    elif family == "inorganic.intermetallic":
        mode = DecompositionMode.METALLIC_SUBLATTICES
        ordered = tuple(sorted(elements))
        for index, first in enumerate(ordered):
            for second in ordered[index + 1:]:
                interactions.append(_interaction((first,), (second,), GrammarOperation.METALLIC_COORDINATION, InteractionPriority.PRIMARY, evidence))
    elif family == "inorganic.oxide":
        mode = DecompositionMode.STRUCTURAL_ANION_SUBSYSTEM
        centres, ligands = tuple(sorted(elements - {"O"})), ("O",)
        interactions.append(_interaction(centres, ligands, GrammarOperation.CENTRE_LIGAND_SHELL, InteractionPriority.PRIMARY, evidence, centres=centres, ligands=ligands))
    elif family == "inorganic.halide":
        mode = DecompositionMode.IONIC_SUBLATTICES
        halogens = elements & reference.chemical.element_set("halogens")
        centres, ligands = tuple(sorted(elements - halogens)), tuple(sorted(halogens))
        interactions.append(_interaction(centres, ligands, GrammarOperation.CENTRE_LIGAND_SHELL, InteractionPriority.PRIMARY, evidence, centres=centres, ligands=ligands))
    elif family in {"inorganic.chalcogenide", "inorganic.nitride", "inorganic.pnictide", "inorganic.carbide", "inorganic.boride", "inorganic.tetrelide"}:
        mode = DecompositionMode.CATION_ANION_SUBSYSTEM
        selectors = {
            "inorganic.chalcogenide": reference.chemical.element_set("chalcogens") - {"O"},
            "inorganic.nitride": frozenset({"N"}),
            "inorganic.pnictide": reference.chemical.element_set("pnictogens") - {"N"},
            "inorganic.carbide": frozenset({"C"}),
            "inorganic.boride": frozenset({"B"}),
            "inorganic.tetrelide": frozenset({"Si", "Ge", "Sn", "Pb"}),
        }
        anions = elements & selectors[family]
        centres, ligands = tuple(sorted(elements - anions)), tuple(sorted(anions))
        interactions.append(_interaction(centres, ligands, GrammarOperation.CENTRE_LIGAND_SHELL, InteractionPriority.PRIMARY, evidence, centres=centres, ligands=ligands))
        if family == "inorganic.chalcogenide":
            interactions.append(_interaction(ligands, ligands, GrammarOperation.INTRA_SUBSYSTEM_BONDS, InteractionPriority.ALLOWED, evidence, centres=ligands, ligands=ligands, layer=InteractionLayer.INTRA_SUBSYSTEM))
    elif family == "organic.molecular":
        mode = DecompositionMode.MOLECULAR_COMPONENTS
        group = tuple(sorted(elements))
        interactions.append(_interaction(group, group, GrammarOperation.COVALENT_NETWORK, InteractionPriority.PRIMARY, evidence, layer=InteractionLayer.INTRAMOLECULAR))
    elif family in {"coordination.metal_organic_discrete", "organic.organometallic_molecular"}:
        mode = DecompositionMode.METAL_ORGANIC
        metals = tuple(sorted(symbol for symbol in elements if reference.elements.by_symbol(symbol).is_metal))
        organic = tuple(sorted(elements - set(metals)))
        interactions.append(_interaction(organic, organic, GrammarOperation.COVALENT_NETWORK, InteractionPriority.PRIMARY, evidence, layer=InteractionLayer.INTRAMOLECULAR))
        donors = tuple(sorted((elements & {"N", "O", "S", "P", "F", "Cl", "Br", "I"}) or (elements & {"C"})))
        interactions.append(_interaction(metals, donors, GrammarOperation.CENTRE_LIGAND_SHELL, InteractionPriority.PRIMARY, evidence, centres=metals, ligands=donors))

    return CompositionGrammar(
        mode=mode,
        candidate_interactions=tuple(interactions),
        confidence=classification.confidence,
        evidence=evidence,
        diagnostics=classification.diagnostics,
        reference_version=reference.chemical.schema_version,
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
