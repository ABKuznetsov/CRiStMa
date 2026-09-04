from __future__ import annotations

import json
import hashlib
from dataclasses import dataclass
from functools import lru_cache
from importlib.resources import files
from types import MappingProxyType
from typing import Any, Mapping


def _freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


@dataclass(frozen=True, slots=True)
class ChemicalReference:
    schema_name: str
    schema_version: str
    resource_sha256: str
    _element_sets: Mapping[str, frozenset[str]]
    _profile_routes: Mapping[str, str]
    _families: Mapping[str, Mapping[str, object]]
    _group_signatures: Mapping[str, Mapping[str, object]]
    _boundary_cases: Mapping[str, Mapping[str, object]]
    _grey_zone_conflicts: tuple[object, ...]
    _structural_contexts: Mapping[str, Mapping[str, object]]
    _grammar_element_sets: Mapping[str, frozenset[str]]
    _composition_family_routes: Mapping[str, object]
    _grammar_templates: Mapping[str, Mapping[str, object]]
    _grammar_routes: Mapping[str, str]

    def element_set(self, identifier: str) -> frozenset[str]:
        return self._element_sets[identifier]

    def profile_route(self, family_identifier: str) -> str:
        return self._profile_routes[family_identifier]

    def family(self, identifier: str) -> Mapping[str, object]:
        return self._families[identifier]

    def group_signature(self, identifier: str) -> Mapping[str, object]:
        return self._group_signatures[identifier]

    def boundary_case(self, identifier: str) -> Mapping[str, object]:
        return self._boundary_cases[identifier]

    def grey_zone_conflicts(self) -> tuple[object, ...]:
        return self._grey_zone_conflicts

    def structural_context(self, identifier: str) -> Mapping[str, object]:
        return self._structural_contexts[identifier]

    def grammar_element_set(self, identifier: str) -> frozenset[str]:
        return self._grammar_element_sets[identifier]

    def composition_family_route(
        self,
        family_identifier: str,
        elements: frozenset[str] | None = None,
    ) -> str:
        route = self._composition_family_routes[family_identifier]
        if isinstance(route, str):
            return route
        if not isinstance(route, Mapping):
            raise ValueError(f"invalid composition-family route: {family_identifier!r}")
        for rule in route.get("rules", ()):
            if not isinstance(rule, Mapping) or rule.get("operator") != "all_other_elements_in":
                raise ValueError(f"unsupported composition routing rule: {rule!r}")
            excluded = frozenset(str(item) for item in rule.get("excluded", ()))
            allowed = self.grammar_element_set(str(rule["set_id"]))
            candidates = (elements or frozenset()) - excluded
            if candidates and candidates <= allowed:
                return str(rule["route"])
        return str(route["default"])

    def grammar_template(self, identifier: str) -> Mapping[str, object]:
        return self._grammar_templates[identifier]

    def grammar_route(self, family_identifier: str) -> str:
        current: str | None = family_identifier
        while current is not None:
            if current in self._grammar_routes:
                return self._grammar_routes[current]
            family = self._families.get(current)
            current = str(family["parent"]) if family is not None and family.get("parent") is not None else None
        raise KeyError(family_identifier)


@dataclass(frozen=True, slots=True)
class ChemicalReferenceIntegrityReport:
    schema_version: str
    resource_sha256: str
    family_count: int
    group_count: int
    boundary_case_count: int
    valid: bool = True


def _assert_acyclic(parents: Mapping[str, str], label: str) -> None:
    for start in parents:
        seen: set[str] = set()
        current: str | None = start
        while current is not None:
            if current in seen:
                raise ValueError(f"cyclic {label} inheritance at {current!r}")
            seen.add(current)
            current = parents.get(current)


def _references(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        return (value,)
    if isinstance(value, tuple) and all(isinstance(item, str) for item in value):
        return value
    raise ValueError("reference field must contain a string or tuple of strings")


def validate_reference_integrity(
    reference: ChemicalReference,
) -> ChemicalReferenceIntegrityReport:
    """Validate semantic references without guessing namespaces from strings."""
    if not reference.schema_version.startswith("3."):
        raise ValueError(f"unsupported chemical reference schema: {reference.schema_version!r}")

    family_ids = set(reference._families)
    group_ids = set(reference._group_signatures)
    route_keys = set(reference._profile_routes)
    if route_keys != family_ids:
        raise ValueError("profile_routing keys must equal family keys")
    registered_profiles = set(reference._profile_routes.values())
    family_parents: dict[str, str] = {}
    for identifier, family in reference._families.items():
        parent = family.get("parent")
        if parent is not None:
            if not isinstance(parent, str) or parent not in family_ids:
                raise ValueError(f"unknown family parent for {identifier!r}: {parent!r}")
            family_parents[identifier] = parent
        profile_id = family.get("profile_id")
        if profile_id != reference._profile_routes[identifier]:
            raise ValueError(f"profile route disagrees with family profile_id for {identifier!r}")
    _assert_acyclic(family_parents, "family")

    operation_ids = {
        "centre_ligand_shell",
        "covalent_network",
        "intra_subsystem_bonds",
        "interstitial_coordination",
        "mixed_anion_coordination",
        "metallic_coordination",
    }
    layer_ids = {
        "structural",
        "interstitial",
        "coordination",
        "intra_subsystem",
        "intramolecular",
        "metallic",
    }
    generic_selectors = {
        "all_elements",
        "remaining_elements",
        "remaining_electropositive_elements",
        "remaining_anion_elements",
        "metal_elements",
        "nonmetal_elements",
        "oxygen",
        "halogens",
    }
    selector_ids = generic_selectors | set(reference._grammar_element_sets)
    for classification_family in reference._composition_family_routes:
        try:
            reference_family = reference.composition_family_route(classification_family)
        except (KeyError, ValueError) as exc:
            raise ValueError(f"invalid composition-family route: {classification_family!r}") from exc
        if not classification_family or reference_family not in family_ids:
            raise ValueError(f"invalid composition-family route: {classification_family!r}")
        route_record = reference._composition_family_routes[classification_family]
        if isinstance(route_record, Mapping):
            for rule in route_record.get("rules", ()):
                target = rule.get("route") if isinstance(rule, Mapping) else None
                set_id = rule.get("set_id") if isinstance(rule, Mapping) else None
                if target not in family_ids or set_id not in reference._grammar_element_sets:
                    raise ValueError(f"invalid composition routing rule for {classification_family!r}")
    for family_id, template_id in reference._grammar_routes.items():
        if family_id not in family_ids:
            raise ValueError(f"unknown grammar-routed family: {family_id!r}")
        if template_id not in reference._grammar_templates:
            raise ValueError(f"unknown grammar template: {template_id!r}")
    for template_id, template in reference._grammar_templates.items():
        interactions = template.get("interactions")
        subsystems = template.get("subsystems")
        if not isinstance(interactions, tuple) or not isinstance(subsystems, tuple):
            raise ValueError(f"grammar template {template_id!r} requires tuple records")
        for record in (*subsystems, *interactions):
            if not isinstance(record, Mapping):
                raise ValueError(f"invalid grammar record in {template_id!r}")
            for field in ("selector", "first", "second", "centres", "ligands"):
                selector = record.get(field)
                if selector is not None and selector not in selector_ids:
                    raise ValueError(f"unknown grammar selector {selector!r} in {template_id!r}")
            operation = record.get("operation")
            if operation is not None and operation not in operation_ids:
                raise ValueError(f"unknown grammar operation {operation!r} in {template_id!r}")
            layer = record.get("layer")
            if layer is not None and layer not in layer_ids:
                raise ValueError(f"unknown grammar layer {layer!r} in {template_id!r}")

    group_parents: dict[str, str] = {}
    for identifier, group in reference._group_signatures.items():
        parent = group.get("parent_group")
        if parent is not None:
            if not isinstance(parent, str) or parent not in group_ids:
                raise ValueError(f"unknown group parent for {identifier!r}: {parent!r}")
            group_parents[identifier] = parent
        family_hint = group.get("family_hint")
        if family_hint is not None and family_hint not in family_ids | registered_profiles:
            raise ValueError(f"unknown family/profile hint for group {identifier!r}: {family_hint!r}")
    _assert_acyclic(group_parents, "group")

    valid_family_or_profile = family_ids | registered_profiles
    for identifier, case in reference._boundary_cases.items():
        composition = case.get("composition")
        refined = case.get("refined")
        variants = case.get("variants")
        for section in (composition, refined):
            if section is None:
                continue
            if not isinstance(section, Mapping):
                raise ValueError(f"boundary section for {identifier!r} must be a mapping")
            for field in ("preferred", "alternatives", "candidates", "preferred_candidates", "competitors", "exclude"):
                for target in _references(section.get(field)):
                    if target not in valid_family_or_profile:
                        raise ValueError(
                            f"unknown family/profile reference in {identifier!r}.{field}: {target!r}"
                        )
            for field in ("recognized_group", "group"):
                target = section.get(field)
                if target is not None and target not in group_ids:
                    raise ValueError(f"unknown group reference in {identifier!r}.{field}: {target!r}")
        if variants is not None:
            if not isinstance(variants, Mapping):
                raise ValueError(f"variants for {identifier!r} must be a mapping")
            for target in variants.values():
                if target not in valid_family_or_profile:
                    raise ValueError(f"unknown variant reference in {identifier!r}: {target!r}")

    return ChemicalReferenceIntegrityReport(
        schema_version=reference.schema_version,
        resource_sha256=reference.resource_sha256,
        family_count=len(family_ids),
        group_count=len(group_ids),
        boundary_case_count=len(reference._boundary_cases),
    )


@lru_cache(maxsize=2)
def load_chemical_reference(version: str = "3.1.0-draft") -> ChemicalReference:
    """Load validated lookup tables; descriptive predicates remain inert data."""
    resource_names = {
        "3.0.0-draft": "chemical_reference_v3.json",
        "3.1.0-draft": "chemical_reference_v3_1.json",
    }
    try:
        resource_name = resource_names[version]
    except KeyError as exc:
        raise ValueError(f"unsupported chemical reference version: {version!r}") from exc
    resource = files("cristma.reference_data").joinpath(f"resources/{resource_name}")
    raw = resource.read_bytes()
    payload = json.loads(raw.decode("utf-8"))
    metadata = payload.get("_meta")
    element_sets = payload.get("element_sets")
    profile_routes = payload.get("profile_routing")
    families = payload.get("families")
    group_signatures = payload.get("group_signatures")
    boundary_cases = payload.get("boundary_cases")
    grey_zone_conflicts = payload.get("grey_zone_conflicts")
    structural_contexts = payload.get("structural_chemistry_contexts")
    grammar_element_sets = payload.get("grammar_element_sets", {})
    composition_family_routes = payload.get("composition_family_routing", {})
    grammar_templates = payload.get("grammar_templates", {})
    grammar_routes = payload.get("grammar_routing", {})
    required_maps = (metadata, element_sets, profile_routes, families, group_signatures, structural_contexts)
    if not all(isinstance(table, dict) for table in required_maps):
        raise ValueError("chemical reference is missing required tables")
    if not isinstance(boundary_cases, list) or not isinstance(grey_zone_conflicts, list):
        raise ValueError("chemical reference is missing required collections")
    schema_version = str(metadata.get("schema_version", ""))
    if not schema_version.startswith("3."):
        raise ValueError(f"unsupported chemical reference schema: {schema_version!r}")
    normalized_sets: dict[str, frozenset[str]] = {}
    for identifier, symbols in element_sets.items():
        if not isinstance(identifier, str) or not isinstance(symbols, list) or not all(
            isinstance(symbol, str) for symbol in symbols
        ):
            raise ValueError("chemical reference contains an invalid element set")
        normalized_sets[identifier] = frozenset(symbols)
    normalized_grammar_sets: dict[str, frozenset[str]] = {}
    for identifier, symbols in grammar_element_sets.items():
        if not isinstance(identifier, str) or not isinstance(symbols, list) or not all(
            isinstance(symbol, str) for symbol in symbols
        ):
            raise ValueError("chemical reference contains an invalid grammar element set")
        normalized_grammar_sets[identifier] = frozenset(symbols)
    if not isinstance(composition_family_routes, dict) or not isinstance(grammar_templates, dict) or not isinstance(grammar_routes, dict):
        raise ValueError("chemical reference contains invalid grammar tables")
    if not all(isinstance(key, str) and isinstance(value, str) for key, value in profile_routes.items()):
        raise ValueError("chemical reference contains invalid profile routing")
    case_index: dict[str, Mapping[str, object]] = {}
    for item in boundary_cases:
        if not isinstance(item, dict) or not isinstance(item.get("id"), str):
            raise ValueError("chemical reference contains an invalid boundary case")
        case_index[item["id"]] = _freeze(item)
    reference = ChemicalReference(
        schema_name=str(metadata.get("schema_name", "")),
        schema_version=schema_version,
        resource_sha256=hashlib.sha256(raw).hexdigest(),
        _element_sets=MappingProxyType(normalized_sets),
        _profile_routes=MappingProxyType(dict(profile_routes)),
        _families=_freeze(families),
        _group_signatures=_freeze(group_signatures),
        _boundary_cases=MappingProxyType(case_index),
        _grey_zone_conflicts=_freeze(grey_zone_conflicts),
        _structural_contexts=_freeze(structural_contexts),
        _grammar_element_sets=MappingProxyType(normalized_grammar_sets),
        _composition_family_routes=MappingProxyType(dict(composition_family_routes)),
        _grammar_templates=_freeze(grammar_templates),
        _grammar_routes=MappingProxyType(dict(grammar_routes)),
    )
    validate_reference_integrity(reference)
    return reference


__all__ = [
    "ChemicalReference",
    "ChemicalReferenceIntegrityReport",
    "load_chemical_reference",
    "validate_reference_integrity",
]
