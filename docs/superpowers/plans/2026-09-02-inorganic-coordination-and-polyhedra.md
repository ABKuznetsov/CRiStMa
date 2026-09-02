# Inorganic Coordination and Polyhedra Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add reproducible inorganic contact resolution, crystallographic-orbit coordination shells, and three-dimensional coordination polyhedra to CRiStMa.

**Architecture:** Convert the existing directed finite/periodic neighbour graphs into canonical `GeometricContact` records, then interpret them with `CompositionGrammar`, covalent radii, and an explicit `ShellResolutionPolicy`. Keep contact selection separate from polyhedron construction: the resolver returns immutable contacts and orbit-projected shells, while `PolyhedronBuilder` accepts only resolved rank-three shells.

**Tech Stack:** Python 3.11+, immutable dataclasses, NumPy 1.26+, pytest 8; no SciPy, spglib, Gemmi, or application dependency at runtime.

**Spec:** `docs/superpowers/specs/2026-09-01-inorganic-coordination-and-polyhedra-design.md`

## Global Constraints

- Preserve `CrystalStructure` and `AtomicView` as canonical immutable inputs; tools return new result objects and never keep a current structure.
- `NeighborFinder` remains chemical-agnostic; all grammar/radius interpretation belongs to `cristma.crystal_chemistry`.
- Use covalent radii for the first normalized-distance implementation: `rho = distance / (r_first + r_second)`.
- Never use expected coordination numbers, compound names, family-specific resolver branches, or a weighted bond score.
- Resolve coordination boundaries collectively per `source_site_id`, then project the decision to each expanded centre.
- Mixed occupancy creates component interpretations, not extra geometric contacts or vertices.
- BVS and geometry evidence are explicit secondary evidence. Until their analyzers exist, return `NOT_AVAILABLE` and `NOT_APPLICABLE`; never simulate evidence.
- A shell is resolvable only when its selected boundary has an observed outer distance group.
- A polyhedron is built only from a `RESOLVED`, affine-rank-three shell.
- Add no runtime dependency beyond existing `numpy>=1.26`.
- Run focused tests after each task. Run the complete suite and wheel smoke test only at the final gate.

## File Map

- `src/cristma/crystallography/local_geometry.py`: canonical finite/periodic geometric contacts and graph conversion.
- `src/cristma/crystal_chemistry/contacts.py`: chemical interpretation, statuses, evidence, shell/result value objects.
- `src/cristma/crystal_chemistry/policy.py`: validated resolution policy and reproducible configuration.
- `src/cristma/crystal_chemistry/resolver.py`: search cutoff, interpretations, grouping, lexicographic boundary selection, orbit projection.
- `src/cristma/crystal_chemistry/polyhedra.py`: rank checks and native small-set convex hull.
- `tests/crystallography/test_local_geometry.py`: contact identity and graph conversion.
- `tests/crystal_chemistry/test_contacts.py`: value-object invariants and mixed occupancy.
- `tests/crystal_chemistry/test_resolver_boundaries.py`: analytic boundary decisions and failure semantics.
- `tests/crystal_chemistry/test_resolver_orbits.py`: orbit consistency and projections.
- `tests/crystal_chemistry/test_polyhedra.py`: hull geometry and non-polyhedral shells.
- `tests/integration/test_inorganic_crystal_chemistry.py`: end-to-end material fixtures.
- `tests/fixtures/crystal_chemistry/`: attributed CIF fixtures and provenance.

---

### Task 1: Canonical geometric contacts

**Files:**
- Create: `src/cristma/crystallography/local_geometry.py`
- Modify: `src/cristma/crystallography/__init__.py`
- Test: `tests/crystallography/test_local_geometry.py`

**Interfaces:**
- Consumes: `AtomicView`, `NeighborGraph`, `PeriodicNeighborGraph`.
- Produces: `GeometricContact` and `geometric_contacts(view, graph) -> tuple[GeometricContact, ...]`.

- [ ] **Step 1: Write failing finite and periodic canonicalization tests**

```python
def test_periodic_reverse_edges_become_one_contact(boundary_view, boundary_graph):
    contacts = geometric_contacts(boundary_view, boundary_graph)
    assert len(contacts) == 1
    assert contacts[0].cell_translation == (-1, 0, 0)
    assert {contacts[0].first_atom_id, contacts[0].second_atom_id} == set(boundary_view.ids)

def test_finite_reverse_edges_become_one_contact(water_view, water_graph):
    contacts = geometric_contacts(water_view, water_graph)
    assert len(contacts) == 2
    assert all(contact.cell_translation is None for contact in contacts)
```

- [ ] **Step 2: Verify the new API is absent**

Run: `pytest tests/crystallography/test_local_geometry.py -v`

Expected: collection fails because `cristma.crystallography.local_geometry` does not exist.

- [ ] **Step 3: Implement the immutable contact and converter**

```python
@dataclass(frozen=True, slots=True)
class GeometricContact:
    contact_id: str
    first_atom_id: str
    second_atom_id: str
    cell_translation: tuple[int, int, int] | None
    distance: float
    vector_cartesian: tuple[float, float, float]
    first_source_site_id: str | None
    second_source_site_id: str | None
    geometric_provenance: str

def geometric_contacts(
    view: AtomicView[AtomicPosition],
    graph: NeighborGraphLike,
) -> tuple[GeometricContact, ...]:
    """Deduplicate directed graph edges using (A, B, t) == (B, A, -t)."""
```

Use the lexicographically smaller forward/reverse key, reverse the vector when the reverse key wins, validate that graph/view atom IDs match, and derive `source_site_id` with `getattr(atom, "source_site_id", None)`.

- [ ] **Step 4: Run focused tests**

Run: `pytest tests/crystallography/test_local_geometry.py tests/geometry/test_finite_neighbors.py tests/geometry/test_periodic_neighbors.py -q`

Expected: all pass; neighbour graph behavior is unchanged.

- [ ] **Step 5: Commit**

```bash
git add src/cristma/crystallography tests/crystallography/test_local_geometry.py
git commit -m "feat: add canonical geometric contacts"
```

### Task 2: Crystal-chemistry contracts and explicit policy

**Files:**
- Create: `src/cristma/crystal_chemistry/__init__.py`
- Create: `src/cristma/crystal_chemistry/contacts.py`
- Create: `src/cristma/crystal_chemistry/policy.py`
- Test: `tests/crystal_chemistry/test_contacts.py`

**Interfaces:**
- Consumes: `GeometricContact`, `GrammarOperation`, `InteractionPriority`, `Diagnostic`.
- Produces: all public immutable result types used by Tasks 3-7.

- [ ] **Step 1: Write failing invariant and configuration tests**

```python
def test_policy_is_explicit_cloneable_and_dimensionless():
    policy = ShellResolutionPolicy(1.45, 0.01, 0.08, 0.01)
    assert policy.get_config() == {
        "candidate_rho_max": 1.45,
        "distance_group_tolerance": 0.01,
        "minimum_shell_gap": 0.08,
        "ambiguity_tolerance": 0.01,
    }
    assert policy.clone(candidate_rho_max=1.60).candidate_rho_max == 1.60
    with pytest.raises(ValueError):
        ShellResolutionPolicy(0.0, 0.01, 0.08, 0.01)

def test_shell_counts_positions_and_occupancy_separately(resolved_contacts):
    shell = CoordinationShell.resolved("site:M", "atom:M", resolved_contacts)
    assert shell.geometric_CN == 4
    assert shell.mean_occupied_neighbors == pytest.approx(3.0)
```

- [ ] **Step 2: Verify tests fail on imports**

Run: `pytest tests/crystal_chemistry/test_contacts.py -v`

Expected: collection fails because `cristma.crystal_chemistry` is absent.

- [ ] **Step 3: Implement the contracts**

```python
class ResolutionStatus(StrEnum):
    RESOLVED = "resolved"
    AMBIGUOUS = "ambiguous"
    INCOMPLETE = "incomplete"
    NOT_APPLICABLE = "not_applicable"

class ContactClassification(StrEnum):
    PRIMARY = "primary"
    SECONDARY = "secondary"

class EvidenceStatus(StrEnum):
    SUPPORTIVE = "supportive"
    NEUTRAL = "neutral"
    CONTRADICTORY = "contradictory"
    NOT_AVAILABLE = "not_available"
    NOT_APPLICABLE = "not_applicable"

@dataclass(frozen=True, slots=True)
class ComponentPairInterpretation:
    first_species: ChemicalSpecies
    second_species: ChemicalSpecies
    first_occupancy: float
    second_occupancy: float
    radius_sum: float
    normalized_distance: float
    occupancy_weight: float
    interaction_type: GrammarOperation
    grammar_priority: InteractionPriority
    centre_elements: tuple[str, ...]
    ligand_elements: tuple[str, ...]

    @property
    def species_symbols(self) -> tuple[str | None, str | None]:
        return (self.first_species.element, self.second_species.element)

@dataclass(frozen=True, slots=True)
class ShellAlternative:
    boundary_group: int
    geometric_CN: int
    relative_gap: float
    internal_spread: float
    strong_contacts_outside: bool

@dataclass(frozen=True, slots=True)
class SecondaryEvidence:
    method: str
    status: EvidenceStatus
    message: str

@dataclass(frozen=True, slots=True)
class ResolvedContact:
    geometric_contact: GeometricContact
    interaction_type: GrammarOperation
    grammar_priority: InteractionPriority
    contact_classification: ContactClassification
    component_interpretations: tuple[ComponentPairInterpretation, ...]
    normalized_distance_min: float
    normalized_distance_max: float
    neighbor_total_occupancy: float
    evidence: tuple[SecondaryEvidence, ...]
    provenance: tuple[tuple[str, object], ...]

@dataclass(frozen=True, slots=True)
class CoordinationShell:
    source_site_id: str
    center_atom_id: str
    contacts: tuple[ResolvedContact, ...]
    geometric_CN: int
    mean_occupied_neighbors: float
    status: ResolutionStatus
    alternatives: tuple[ShellAlternative, ...]
    evidence: tuple[SecondaryEvidence, ...]
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]

    @property
    def diagnostic_codes(self) -> tuple[str, ...]:
        return tuple(item.code for item in self.diagnostics)

    @classmethod
    def resolved(
        cls,
        source_site_id: str,
        center_atom_id: str,
        contacts: tuple[ResolvedContact, ...],
    ) -> CoordinationShell:
        return cls(
            source_site_id=source_site_id,
            center_atom_id=center_atom_id,
            contacts=contacts,
            geometric_CN=len(contacts),
            mean_occupied_neighbors=math.fsum(
                item.neighbor_total_occupancy for item in contacts
            ),
            status=ResolutionStatus.RESOLVED,
            alternatives=(), evidence=(), diagnostics=(), provenance=(),
        )

@dataclass(frozen=True, slots=True)
class CrystalChemistryResolution:
    contacts: tuple[ResolvedContact, ...]
    coordination_shells: tuple[CoordinationShell, ...]
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]

    @property
    def diagnostic_codes(self) -> tuple[str, ...]:
        return tuple(item.code for item in self.diagnostics)
```

Store provenance as immutable tuples of `(key, value)` pairs. `CoordinationShell` validates that computed counts match its contacts; do not accept independently supplied contradictory counts.

- [ ] **Step 4: Implement `ShellResolutionPolicy`**

```python
@dataclass(frozen=True, slots=True)
class ShellResolutionPolicy:
    candidate_rho_max: float
    distance_group_tolerance: float
    minimum_shell_gap: float
    ambiguity_tolerance: float

    def get_config(self) -> dict[str, float]:
        return {
            "candidate_rho_max": self.candidate_rho_max,
            "distance_group_tolerance": self.distance_group_tolerance,
            "minimum_shell_gap": self.minimum_shell_gap,
            "ambiguity_tolerance": self.ambiguity_tolerance,
        }

    def clone(self, **changes: float) -> ShellResolutionPolicy:
        return replace(self, **changes)
```

Reject booleans, non-finite values, and values `<= 0`. Do not add a module-level default.

- [ ] **Step 5: Run tests and commit**

Run: `pytest tests/crystal_chemistry/test_contacts.py -q`

```bash
git add src/cristma/crystal_chemistry tests/crystal_chemistry/test_contacts.py
git commit -m "feat: define crystal chemistry result contracts"
```

### Task 3: Grammar coverage and component interpretations

**Files:**
- Modify: `src/cristma/chemistry/grammar.py`
- Create: `src/cristma/crystal_chemistry/resolver.py`
- Modify: `tests/chemistry/test_grammar.py`
- Test: `tests/crystal_chemistry/test_interpretations.py`

**Interfaces:**
- Consumes: `CompositionGrammar`, `ReferenceData.covalent_radii`, `GeometricContact`, endpoint components.
- Produces: `derive_search_cutoff(...)` and private `_interpret_contact(...)` used by the resolver.

- [ ] **Step 1: Add grammar tests for all directed operation names**

```python
def test_directed_coordination_operations_are_distinct():
    assert GrammarOperation.INTERSTITIAL_COORDINATION.value == "interstitial_coordination"
    assert GrammarOperation.MIXED_ANION_COORDINATION.value == "mixed_anion_coordination"
```

Extend `GrammarOperation` with those values without changing existing family compilation yet.

- [ ] **Step 2: Write interpretation tests**

```python
def test_mixed_sites_make_four_interpretations_not_four_contacts(mixed_contact_input):
    outcome = _interpret_contact(*mixed_contact_input)
    assert len(outcome.interpretations) == 4
    assert {pair.species_symbols for pair in outcome.interpretations} == {
        ("Ca", "O"), ("Ca", "F"), ("Sr", "O"), ("Sr", "F")
    }

def test_search_cutoff_uses_largest_allowed_radius_sum(grammar, reference, policy):
    assert derive_search_cutoff(grammar, reference, policy) == pytest.approx(
        max_allowed_radius_sum(grammar, reference) * policy.candidate_rho_max
    )
```

- [ ] **Step 3: Run tests to observe missing functions**

Run: `pytest tests/chemistry/test_grammar.py tests/crystal_chemistry/test_interpretations.py -q`

Expected: grammar enum test and resolver imports fail.

- [ ] **Step 4: Implement exact pair matching and normalization**

Implement helpers with these signatures:

```python
def derive_search_cutoff(
    grammar: CompositionGrammar,
    reference: ReferenceData,
    policy: ShellResolutionPolicy,
) -> float:
    sums = _allowed_radius_sums(grammar, reference)
    if not sums:
        raise ValueError("grammar has no component pairs with known covalent radii")
    return max(sums) * policy.candidate_rho_max

def _interpret_contact(
    contact: GeometricContact,
    first_components: tuple[SiteComponent, ...],
    second_components: tuple[SiteComponent, ...],
    grammar: CompositionGrammar,
    reference: ReferenceData,
    policy: ShellResolutionPolicy,
) -> InterpretationOutcome:
    records: list[ComponentPairInterpretation] = []
    for first in first_components:
        for second in second_components:
            requests = _matching_interactions(first, second, grammar)
            if not requests:
                continue
            radius_sum = (
                reference.covalent_radii.find(first.element).value
                + reference.covalent_radii.find(second.element).value
            )
            rho = contact.distance / radius_sum
            if rho <= policy.candidate_rho_max:
                for request in requests:
                    records.append(ComponentPairInterpretation(
                        first.species, second.species,
                        float(first.occupancy.value), float(second.occupancy.value),
                        radius_sum, rho,
                        float(first.occupancy.value) * float(second.occupancy.value),
                        request.operation, request.priority,
                        request.centre_elements, request.ligand_elements,
                    ))
    return InterpretationOutcome(tuple(records), ())
```

Pair matching is unordered for network contacts but respects `centre_elements`/`ligand_elements` when shells are projected. For each permitted pair, look up both radii, compute `radius_sum`, `rho`, and the occupancy product. Exclude `rho > candidate_rho_max`. Return machine diagnostics for missing radii/unknown species through an internal interpretation outcome rather than raising.

Define the internal result beside the helpers, then implement the helpers used above with exact signatures:

```python
@dataclass(frozen=True, slots=True)
class InterpretationOutcome:
    interpretations: tuple[ComponentPairInterpretation, ...]
    diagnostics: tuple[Diagnostic, ...]

def _allowed_radius_sums(
    grammar: CompositionGrammar,
    reference: ReferenceData,
) -> tuple[float, ...]:
    values: list[float] = []
    for request in grammar.candidate_interactions:
        for first in request.first_elements:
            for second in request.second_elements:
                try:
                    values.append(
                        reference.covalent_radii.find(first).value
                        + reference.covalent_radii.find(second).value
                    )
                except KeyError:
                    continue
    return tuple(values)

def _matching_interactions(
    first: SiteComponent,
    second: SiteComponent,
    grammar: CompositionGrammar,
) -> tuple[CandidateInteraction, ...]:
    first_symbol, second_symbol = first.element, second.element
    if first_symbol is None or second_symbol is None:
        return ()
    matches: list[CandidateInteraction] = []
    for request in grammar.candidate_interactions:
        forward = (
            first_symbol in request.first_elements
            and second_symbol in request.second_elements
        )
        reverse = (
            second_symbol in request.first_elements
            and first_symbol in request.second_elements
        )
        if forward or reverse:
            matches.append(request)
    return tuple(matches)
```

Their bodies perform only catalog lookup and exact element-set matching; they do not inspect coordinates, compound names, or expected CN. Each match remains a separate scope keyed by source site, interaction, and centre view. Missing radii may be omitted only while deriving the numerical cutoff: interpretation records `radius_missing`, and a missing PRIMARY pair makes its scope `INCOMPLETE`.

- [ ] **Step 5: Add explicit reference/source failure tests**

```python
@pytest.mark.parametrize("case,code", [
    ("missing_radius", "crystal_chemistry.contact.radius_missing"),
    ("unknown_species", "crystal_chemistry.contact.species_unknown"),
    ("invalid_source_occupancy", "crystal_chemistry.contact.occupancy_invalid"),
])
def test_unusable_component_data_is_diagnostic_not_a_guessed_contact(case, code):
    outcome = interpret_source_case(case)
    assert not outcome.interpretations
    assert code in {item.code for item in outcome.diagnostics}
```

The invalid-occupancy fixture is an internal mapper outcome, because canonical `SiteComponent` correctly rejects invalid occupancy before resolver entry.

- [ ] **Step 6: Run tests and commit**

Run: `pytest tests/chemistry/test_grammar.py tests/crystal_chemistry/test_interpretations.py tests/reference_data/test_radii.py -q`

```bash
git add src/cristma/chemistry/grammar.py src/cristma/crystal_chemistry/resolver.py tests/chemistry/test_grammar.py tests/crystal_chemistry/test_interpretations.py
git commit -m "feat: interpret inorganic contact components"
```

### Task 4: Traceable shell-boundary selection

**Files:**
- Modify: `src/cristma/crystal_chemistry/resolver.py`
- Test: `tests/crystal_chemistry/test_resolver_boundaries.py`

**Interfaces:**
- Consumes: normalized contact interpretations and `ShellResolutionPolicy`.
- Produces: private `_group_distances`, `_candidate_boundaries`, `_select_boundary` and public evidence in `ShellAlternative`.

- [ ] **Step 1: Write analytic tests for grouping and formulas**

```python
def test_candidate_uses_relative_gap_and_internal_spread(policy):
    decision = resolve_rho_fixture((0.98, 1.00, 1.01, 1.28), policy)
    assert decision.status is ResolutionStatus.RESOLVED
    assert decision.selected.geometric_CN == 3
    assert decision.selected.relative_gap == pytest.approx((1.28 - 1.01) / 1.01)
    assert decision.selected.internal_spread == pytest.approx((1.01 - 0.98) / 1.00)

def test_close_non_dominated_boundaries_remain_ambiguous(policy):
    decision = resolve_rho_fixture((0.95, 1.00, 1.10, 1.16, 1.28), policy)
    assert decision.status is ResolutionStatus.AMBIGUOUS
    assert len(decision.alternatives) >= 2
```

- [ ] **Step 2: Add completeness and strong-outside tests**

```python
def test_last_group_cannot_prove_a_shell_boundary(policy):
    decision = resolve_rho_fixture((0.98, 1.00, 1.01), policy)
    assert decision.status is ResolutionStatus.INCOMPLETE
    assert diagnostic_codes(decision) == {
        "crystal_chemistry.shell.search_boundary_not_observed"
    }

def test_too_few_groups_is_incomplete(policy):
    decision = resolve_rho_fixture((1.00,), policy)
    assert "crystal_chemistry.shell.candidates_insufficient" in diagnostic_codes(decision)

def test_non_dominated_boundaries_report_stable_code(policy):
    decision = resolve_rho_fixture((0.95, 1.00, 1.10, 1.16, 1.28), policy)
    assert "crystal_chemistry.shell.boundary_ambiguous" in diagnostic_codes(decision)
```

- [ ] **Step 3: Verify failures**

Run: `pytest tests/crystal_chemistry/test_resolver_boundaries.py -q`

- [ ] **Step 4: Implement hard filters, grouping and lexicographic comparison**

Use sorted normalized-distance groups, grouping adjacent values whose dimensionless difference is within `distance_group_tolerance`. For every boundary with an outer group, calculate the exact spec formulas. Compare gap, spread, strong-outside, then grammar priority; apply `ambiguity_tolerance` separately at each numeric comparison. Return all non-dominated alternatives. Add fixed BVS/geometry evidence records with `NOT_AVAILABLE`/`NOT_APPLICABLE` only.

- [ ] **Step 5: Run and commit**

Run: `pytest tests/crystal_chemistry/test_resolver_boundaries.py -q`

```bash
git add src/cristma/crystal_chemistry/resolver.py tests/crystal_chemistry/test_resolver_boundaries.py
git commit -m "feat: resolve normalized shell boundaries"
```

### Task 5: Public resolver, orbit consistency, and shell projection

**Files:**
- Modify: `src/cristma/crystal_chemistry/resolver.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Test: `tests/crystal_chemistry/test_resolver_orbits.py`

**Interfaces:**
- Consumes: `CrystalStructure`, `CompositionGrammar`, explicit policy and reference data.
- Produces: `CoordinationShellResolver.resolve(...) -> CrystalChemistryResolution`.

- [ ] **Step 1: Write end-to-end resolver API test**

```python
def test_resolver_derives_cutoff_and_returns_reproducible_provenance(naf_crystal, naf_grammar):
    result = CoordinationShellResolver(policy=POLICY).resolve(naf_crystal, naf_grammar)
    provenance = dict(result.provenance)
    assert provenance["policy"] == POLICY.get_config()
    assert provenance["search_cutoff_angstrom"] > 0
    assert provenance["grammar_method"] == "cristma.composition_grammar:1"
    assert provenance["resolver_method"] == "cristma.coordination_shell_resolver:1"
```

- [ ] **Step 2: Write orbit projection and inconsistency tests**

```python
def test_one_orbit_decision_projects_to_every_equivalent_center(equivalent_centers):
    result = resolve(equivalent_centers)
    shells = [s for s in result.coordination_shells if s.source_site_id == "site:M"]
    expected = [
        atom for atom in equivalent_centers.atomic_view().atoms
        if atom.source_site_id == "site:M"
    ]
    assert len(shells) == len(expected)
    assert len({shell.geometric_CN for shell in shells}) == 1

def test_orbit_signature_mismatch_is_incomplete(inconsistent_graph):
    result = resolve_with_graph(inconsistent_graph)
    assert all(shell.status is ResolutionStatus.INCOMPLETE for shell in result.coordination_shells)
    assert "crystal_chemistry.shell.symmetry_inconsistent" in result.diagnostic_codes
```

- [ ] **Step 3: Verify failures**

Run: `pytest tests/crystal_chemistry/test_resolver_orbits.py -q`

- [ ] **Step 4: Implement resolver orchestration**

```python
@dataclass(frozen=True, slots=True)
class CoordinationShellResolver:
    policy: ShellResolutionPolicy
    reference: ReferenceData = field(default_factory=ReferenceData.default)

    def resolve(
        self,
        structure: CrystalStructure,
        grammar: CompositionGrammar,
    ) -> CrystalChemistryResolution:
        view = structure.atomic_view()
        cutoff = derive_search_cutoff(grammar, self.reference, self.policy)
        graph = NeighborFinder(cutoff=cutoff).find(view)
        contacts = geometric_contacts(view, graph)
        interpreted = self._resolve_contacts(view, contacts, grammar)
        shells = self._resolve_orbits(view, interpreted, grammar)
        return self._make_result(structure, grammar, cutoff, interpreted, shells)
```

The method expands the structure, derives cutoff, calls `NeighborFinder`, canonicalizes contacts, interprets components, resolves directed shell operations collectively by source site, and keeps network operations only in `result.contacts`. Orbit signatures use neighbour source site, component pair, interaction, rho group, and multiplicity—never raw expanded IDs or translation vectors.

- [ ] **Step 5: Add mixed-occupancy disagreement and vacancy tests**

```python
def test_component_boundaries_that_disagree_are_ambiguous(mixed_boundary_crystal):
    shell = only_shell(resolve(mixed_boundary_crystal))
    assert shell.status is ResolutionStatus.AMBIGUOUS
    assert "crystal_chemistry.shell.mixed_occupancy_disagreement" in shell.diagnostic_codes

def test_vacancy_changes_mean_occupancy_not_geometric_cn(vacancy_crystal):
    shell = only_shell(resolve(vacancy_crystal))
    assert shell.geometric_CN == 4
    assert shell.mean_occupied_neighbors == pytest.approx(3.0)
```

- [ ] **Step 6: Run focused resolver suite and commit**

Run: `pytest tests/crystal_chemistry/test_contacts.py tests/crystal_chemistry/test_interpretations.py tests/crystal_chemistry/test_resolver_boundaries.py tests/crystal_chemistry/test_resolver_orbits.py -q`

```bash
git add src/cristma/crystal_chemistry tests/crystal_chemistry
git commit -m "feat: resolve orbit-level coordination shells"
```

### Task 6: Native three-dimensional polyhedron construction

**Files:**
- Create: `src/cristma/crystal_chemistry/polyhedra.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Test: `tests/crystal_chemistry/test_polyhedra.py`

**Interfaces:**
- Consumes: one `CoordinationShell` and the source `AtomicView`.
- Produces: `CoordinationPolyhedron`, `PolyhedronBuildResult`, `PolyhedronBuilder.build(...)`.

- [ ] **Step 1: Write tetrahedron and octahedron tests**

```python
def test_tetrahedral_shell_builds_closed_hull(tetrahedral_shell, view):
    result = PolyhedronBuilder().build(tetrahedral_shell, view)
    assert result.status is ResolutionStatus.RESOLVED
    assert len(result.polyhedron.faces) == 4
    assert result.polyhedron.volume > 0
    assert result.polyhedron.center_offset == pytest.approx(0.0)

def test_octahedral_shell_has_eight_triangular_faces(octahedral_shell, view):
    polyhedron = PolyhedronBuilder().build(octahedral_shell, view).polyhedron
    assert len(polyhedron.faces) == 8
```

- [ ] **Step 2: Write rejection tests**

```python
@pytest.mark.parametrize("fixture,status", [
    ("linear", ResolutionStatus.NOT_APPLICABLE),
    ("planar", ResolutionStatus.NOT_APPLICABLE),
    ("ambiguous", ResolutionStatus.AMBIGUOUS),
    ("incomplete", ResolutionStatus.INCOMPLETE),
])
def test_non_success_reason_is_preserved(fixture, status, request):
    shell, view = request.getfixturevalue(fixture)
    result = PolyhedronBuilder().build(shell, view)
    assert result.polyhedron is None
    assert result.status is status
```

- [ ] **Step 3: Verify failures**

Run: `pytest tests/crystal_chemistry/test_polyhedra.py -q`

- [ ] **Step 4: Implement a native small-set hull**

```python
@dataclass(frozen=True, slots=True)
class CoordinationPolyhedron:
    polyhedron_id: str
    source_site_id: str
    center_atom_id: str
    shell_provenance: tuple[tuple[str, object], ...]
    vertex_contacts: tuple[ResolvedContact, ...]
    local_vertices: tuple[tuple[float, float, float], ...]
    faces: tuple[tuple[int, ...], ...]
    volume: float
    geometric_centroid: tuple[float, float, float]
    center_offset: float
    diagnostics: tuple[Diagnostic, ...] = ()

@dataclass(frozen=True, slots=True)
class PolyhedronBuildResult:
    status: ResolutionStatus
    polyhedron: CoordinationPolyhedron | None
    diagnostics: tuple[Diagnostic, ...]

class PolyhedronBuilder:
    def build(
        self,
        shell: CoordinationShell,
        view: AtomicView,
    ) -> PolyhedronBuildResult:
        if shell.status is not ResolutionStatus.RESOLVED:
            return PolyhedronBuildResult.from_shell_failure(shell)
        vertices = _local_vertices(shell, view)
        if np.linalg.matrix_rank(vertices[1:] - vertices[0]) < 3:
            return PolyhedronBuildResult.not_applicable(
                "crystal_chemistry.polyhedron.not_three_dimensional"
            )
        faces = _convex_hull_faces(vertices)
        return PolyhedronBuildResult.resolved(
            _make_polyhedron(shell, vertices, faces)
        )

    def validate_orbit(
        self,
        polyhedra: tuple[CoordinationPolyhedron, ...],
    ) -> tuple[Diagnostic, ...]:
        signatures = {_face_topology_signature(item) for item in polyhedra}
        if len(signatures) <= 1:
            return ()
        return (Diagnostic(
            Severity.ERROR,
            "crystal_chemistry.polyhedron.symmetry_inconsistent",
            "symmetry-equivalent centres have different face topology",
        ),)
```

Implement `PolyhedronBuildResult.resolved(...)` and `.not_applicable(code)` as constructors that guarantee exactly one of a polyhedron or failure diagnostics. Implement `_local_vertices`, `_convex_hull_faces`, and `_make_polyhedron` in the same module; they are private numerical helpers, not public APIs.

Translate periodic ligand positions using each contact vector. Check affine rank with `numpy.linalg.matrix_rank`. Enumerate all vertex triples; retain supporting planes with all remaining points on one side; merge coplanar triangles into maximal ordered polygon faces; orient them outward; calculate volume and centroid using an internal triangulation. Keep one vertex per geometric contact regardless of occupancy.

Add a cube test: eight vertices expose six quadrilateral faces, not twelve triangles.

- [ ] **Step 5: Add symmetry-topology test, run and commit**

```python
def test_equivalent_centres_have_same_orbit_aware_face_topology(equivalent_polyhedra):
    signatures = {polyhedron_face_signature(item) for item in equivalent_polyhedra}
    assert len(signatures) == 1

def test_inconsistent_equivalent_topology_returns_diagnostic(inconsistent_polyhedra):
    diagnostics = PolyhedronBuilder().validate_orbit(inconsistent_polyhedra)
    assert "crystal_chemistry.polyhedron.symmetry_inconsistent" in {
        item.code for item in diagnostics
    }
```

Run: `pytest tests/crystal_chemistry/test_polyhedra.py -q`

```bash
git add src/cristma/crystal_chemistry/polyhedra.py src/cristma/crystal_chemistry/__init__.py tests/crystal_chemistry/test_polyhedra.py
git commit -m "feat: build coordination polyhedra"
```

### Task 7: Acceptance fixture corpus and material-independent integration

**Files:**
- Create: `tests/fixtures/crystal_chemistry/PROVENANCE.md`
- Create: attributed CIF files under `tests/fixtures/crystal_chemistry/`
- Create: `tests/integration/test_inorganic_crystal_chemistry.py`

**Interfaces:**
- Consumes: public `cristma.read`, `ChemistryAnalyzer`, `CoordinationShellResolver`, `PolyhedronBuilder`.
- Produces: calibrated acceptance evidence without embedding expected CN in production code.

- [ ] **Step 1: Commit fixture provenance before assertions**

For every NaF, SiC, Si3N4, FeS2, CaN2, Na3P, Bi2Te3, CaMoO4, LiB3O5, and anorthite fixture, record source URL/database identifier, license/redistribution basis, retrieval date, formula, space group, and checksum in `PROVENANCE.md`. Prefer already committed attributed CIFs where chemically suitable; do not copy unattributed structures.

- [ ] **Step 2: Add network/contact acceptance tests**

```python
@pytest.mark.parametrize("name,operation", [
    ("SiC", GrammarOperation.COVALENT_NETWORK),
    ("FeS2", GrammarOperation.INTRA_SUBSYSTEM_BONDS),
    ("NiAl", GrammarOperation.METALLIC_COORDINATION),
])
def test_network_materials_return_contacts_without_forced_polyhedra(name, operation):
    result = calculate(name)
    assert any(contact.interaction_type is operation for contact in result.contacts)
```

If no redistributable NiAl fixture is obtained, cover metallic contacts with an attributed analytic B2 fixture documented as generated coordinates.

- [ ] **Step 3: Add shell/polyhedron acceptance tests**

```python
@pytest.mark.parametrize("name,center,cn", [
    ("NaF", "Na", 6), ("Si3N4", "Si", 4), ("FeS2", "Fe", 6),
    ("CaMoO4", "Mo", 4),
])
def test_established_shells_are_discovered_not_supplied(name, center, cn):
    result = calculate(name)
    shells = resolved_shells_for_element(result, center)
    assert shells and {shell.geometric_CN for shell in shells} == {cn}
```

Add explicit tests for Ca-N plus N-N where present, Na-P shells, retained Bi-Te secondary candidates, and BO3/BO4 in LiB3O5. For anorthite accept resolved Al/SiO4 polyhedra or an explicit `AMBIGUOUS`/`INCOMPLETE` shell with no forced polyhedron. The helper passes only structure, grammar, reference, and policy to production code.

- [ ] **Step 4: Calibrate one explicit test policy, not a hidden default**

Keep `ACCEPTANCE_POLICY = ShellResolutionPolicy(1.60, 0.01, 0.08, 0.01)` in the test module as the first calibration candidate. If one policy cannot cover the corpus, report failing structures and revise the scientific algorithm; do not add compound-name branches or per-compound production constants.

- [ ] **Step 5: Run only the acceptance slice and commit**

Run: `pytest tests/integration/test_inorganic_crystal_chemistry.py -q`

```bash
git add tests/fixtures/crystal_chemistry tests/integration/test_inorganic_crystal_chemistry.py
git commit -m "test: add inorganic crystal chemistry corpus"
```

### Task 8: Public package surface, documentation, and final verification

**Files:**
- Modify: `src/cristma/crystallography/__init__.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Modify: `README.md`
- Create: `docs/inorganic-crystal-chemistry.md`
- Test: `tests/integration/test_public_api.py`

**Interfaces:**
- Consumes: all completed scientific types/tools.
- Produces: stable internal-development import surface and runnable examples.

- [ ] **Step 1: Add public import smoke test**

```python
def test_crystal_chemistry_tools_import_from_package():
    from cristma.crystallography import GeometricContact, geometric_contacts
    from cristma.crystal_chemistry import (
        CoordinationShellResolver, PolyhedronBuilder, ShellResolutionPolicy,
    )
```

- [ ] **Step 2: Document the minimal workflow and scientific limits**

```python
structure = cristma.read("sample.cif").structures[0]
chemistry = ChemistryAnalyzer().analyze(Composition.from_structure(structure))
policy = ShellResolutionPolicy(1.45, 0.01, 0.08, 0.01)
resolution = CoordinationShellResolver(policy).resolve(structure, chemistry.grammar)
polyhedra = [PolyhedronBuilder().build(shell, structure.atomic_view())
             for shell in resolution.coordination_shells]
```

State that the numbers are explicit example configuration, not a universal preset; document status meanings, missing BVS/geometry analyzers, and provenance fields.

- [ ] **Step 3: Run focused package tests**

Run: `pytest tests/crystallography/test_local_geometry.py tests/crystal_chemistry tests/integration/test_inorganic_crystal_chemistry.py tests/integration/test_public_api.py -q`

Expected: all pass.

- [ ] **Step 4: Run the complete suite once**

Run: `pytest -q`

Expected: all existing and new tests pass.

- [ ] **Step 5: Build and inspect a wheel in a temporary directory**

Run: `python -m build --wheel --outdir /tmp/cristma-dist`

Run: `python -c "import zipfile, pathlib; p=next(pathlib.Path('/tmp/cristma-dist').glob('*.whl')); z=zipfile.ZipFile(p); assert any(n.endswith('crystal_chemistry/resolver.py') for n in z.namelist()); print(p)"`

Expected: wheel contains the new package and package metadata still declares only NumPy as a base runtime dependency.

- [ ] **Step 6: Commit final surface**

```bash
git add src/cristma README.md docs/inorganic-crystal-chemistry.md tests/integration/test_public_api.py
git commit -m "docs: expose inorganic crystal chemistry workflow"
```

## Final Review Checklist

- Every stable diagnostic from the spec has a direct test.
- No production source contains material names from the acceptance corpus.
- No expected coordination-number table or weighted score exists.
- Mixed occupancy keeps one edge/vertex and reports both geometric CN and mean occupied neighbours.
- Every resolved boundary has an observed outer group.
- Orbit inconsistency returns `INCOMPLETE`, not `AMBIGUOUS`.
- Network operations remain contacts; only directed coordination operations produce shells.
- Polyhedra accept only resolved rank-three shells.
- BVS/geometry evidence is explicit and honestly unavailable/not applicable.
- Provenance records policy, cutoff, maximum observed rho, outer-group completeness, grammar/reference/resolver versions, and structure identity.
- Full suite and built-wheel smoke test pass.
