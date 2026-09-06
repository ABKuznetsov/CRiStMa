# Orbit-First Crystal Chemistry Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace expanded-first crystal chemistry with validated asymmetric-unit pair orbits, orbit-first chemical and shell resolution, and outward-only materialization.

**Architecture:** Direct-space analysis starts from independent sites and a validated `SymmetryContext`. Exact affine-periodic relations and endpoint stabilizers produce canonical geometric pair orbits; chemistry, oriented incidences, shells, polyhedra, and structural hierarchy consume those orbits. Expanded contacts are derived only by CrIStMa materializers, and the legacy expanded-first pipeline is deleted at cutover.

**Tech Stack:** Python 3.11+, frozen dataclasses with slots, `fractions.Fraction`, NumPy, SciPy only where existing polyhedron geometry already requires it, pytest, native CrIStMa CIF/structure/symmetry models.

**Spec:** `docs/superpowers/specs/2026-09-06-orbit-first-crystal-chemistry-design.md`

## Global Constraints

- CrIStMa contains scientific calculations only; no CRAFT, Finder, scene, colour, visibility, table, matching, or refinement logic.
- Direct-space orbit-first analysis accepts only a validated `SymmetryContext`.
- Diffraction v1 continues to require a catalog `SpaceGroupSetting` and is not changed by this plan.
- Exact symmetry and lattice-relation algebra uses `Fraction`; numerical tolerance is restricted to coordinates, cell-metric compatibility, and cutoff boundaries.
- Operation, relation, and orbit identities are deterministic and independent of CIF operation order.
- `geometry_orbit_id` excludes cutoff, distance, and metric cell parameters.
- Pair multiplicity and local incidence multiplicity remain distinct.
- `PRIMARY` and `SECONDARY` belong to shell alternatives, never to incidence identity.
- Materialized contacts never feed back into scientific analysis.
- No new runtime dependency is added.
- The new public API is not released until downstream hierarchy consumers are migrated and the legacy expanded-first route is removed.
- Tests and development specifications are not package data and therefore do not enter wheels or sdists unless packaging configuration is explicitly changed.

## File Structure

New focused modules:

- `src/cristma/crystallography/symmetry_context.py` — context construction, validation, canonical operation identity, and context fingerprints.
- `src/cristma/crystallography/periodic_relation.py` — exact affine-periodic relation algebra.
- `src/cristma/crystallography/asu_mapping.py` — independent-site stabilizers, site images, and mapping fingerprints.
- `src/cristma/crystallography/symmetry_pairs.py` — cutoff-complete buffered pair search and geometric pair-table models.
- `src/cristma/crystallography/pair_canonical.py` — stabilizer double-coset canonicalization, endpoint exchange, ownership, and multiplicity.
- `src/cristma/crystal_chemistry/orbit_contacts.py` — chemical interpretations of geometric pair orbits.
- `src/cristma/crystal_chemistry/incidence_orbits.py` — oriented local incidences and effective participating occupancy.
- `src/cristma/crystal_chemistry/shell_orbits.py` — shell alternatives, shell resolution, and aggregate status.
- `src/cristma/crystal_chemistry/contact_analysis.py` — orchestration and primary `ContactAnalysisResult`.
- `src/cristma/crystal_chemistry/materialize.py` — reference-cell/cell-range materialization and compatibility views.

Existing modules modified during migration:

- `src/cristma/crystallography/__init__.py` — intentional public exports for symmetry context and pair geometry.
- `src/cristma/crystal_chemistry/__init__.py` — intentional public exports for orbit-first chemistry.
- `src/cristma/crystal_chemistry/polyhedra.py` — representative local geometry from shell orbits.
- `src/cristma/crystal_chemistry/polyhedron_orbits.py` — orbit identity from orbit-first shell inputs.
- `src/cristma/crystal_chemistry/structural_units.py` — units consume polyhedron/shell orbits.
- `src/cristma/crystal_chemistry/structural_graph.py` — graph edges reference resolved orbit identities.
- `src/cristma/crystal_chemistry/representation.py` — representation selection operates on orbit graph records.
- `src/cristma/crystal_chemistry/periodic_connectivity.py` — translation rank from orbit relations.
- `src/cristma/crystal_chemistry/structural_blocks.py` — blocks and block orbits from orbit-first components.
- `src/cristma/crystal_chemistry/ring_finder.py`, `_ring_search.py`, `_ring_symmetry.py`, and `rings.py` — ring search and identity from orbit relations without expanded contacts.
- `tools/smoke_cif_corpus.py` — acceptance route uses `ContactAnalyzer` and explicit materialization only for reported instance counts.

Legacy modules removed at cutover:

- `src/cristma/crystal_chemistry/resolver.py`;
- `src/cristma/crystal_chemistry/contact_orbits.py`;
- expanded-first models in `src/cristma/crystal_chemistry/contacts.py` after their orbit-first replacements are exported.

---

### Task 1: Validated `SymmetryContext`

**Files:**
- Create: `src/cristma/crystallography/symmetry_context.py`
- Modify: `src/cristma/crystallography/__init__.py`
- Test: `tests/orbit_first/test_symmetry_context.py`

**Interfaces:**
- Consumes: `UnitCell`, `AffineOperation`, `SpaceGroupDefinition`, and `SpaceGroupSetting`.
- Produces: `SymmetryContext`, `SymmetrySourceKind`, `DirectBasisConvention`, `SymmetryContextInvariantError`, `canonical_operation_key(operation) -> str`.

- [ ] **Step 1: Write failing context-construction and validation tests**

```python
def test_explicit_operation_order_does_not_change_context_identity():
    first = SymmetryContext.from_operations((identity, inversion), cell)
    second = SymmetryContext.from_operations((inversion, identity), cell)
    assert first.operation_keys == second.operation_keys
    assert first.symmetry_action_fingerprint == second.symmetry_action_fingerprint

def test_valid_explicit_operations_are_complete_without_setting_id():
    context = SymmetryContext.from_operations((identity, inversion), cell)
    assert context.source_kind is SymmetrySourceKind.VALID_EXPLICIT_OPERATIONS
    assert context.setting_id is None

def test_non_closed_operation_set_is_rejected():
    with pytest.raises(SymmetryContextInvariantError) as caught:
        SymmetryContext.from_operations((identity, quarter_turn), square_cell)
    assert caught.value.code == "symmetry.context.group_invalid"
```

- [ ] **Step 2: Run the tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_symmetry_context.py`

Expected: collection fails because `symmetry_context` and its public types do not exist.

- [ ] **Step 3: Implement exact canonicalization and validation**

```python
class SymmetrySourceKind(StrEnum):
    CATALOG_SETTING = "catalog_setting"
    VALID_EXPLICIT_OPERATIONS = "valid_explicit_operations"
    EXPLICIT_IDENTITY_FALLBACK = "explicit_identity_fallback"

@dataclass(frozen=True, slots=True)
class SymmetryContext:
    operations: tuple[AffineOperation, ...]
    operation_keys: tuple[str, ...]
    basis_convention: DirectBasisConvention
    cell_fingerprint: str
    symmetry_action_fingerprint: str
    setting_id: str | None
    source_kind: SymmetrySourceKind
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]

    @classmethod
    def from_operations(cls, operations, cell, *, provenance=()):
        canonical = validate_and_canonicalize_operations(operations, cell)
        return build_context(
            canonical, cell,
            source_kind=SymmetrySourceKind.VALID_EXPLICIT_OPERATIONS,
            setting_id=None,
            provenance=provenance,
        )

    @classmethod
    def from_definition(cls, definition, cell):
        return cls.from_operations(
            definition.operations,
            cell,
            provenance=(("source", "space_group_definition"),),
        )

    @classmethod
    def from_setting(cls, setting, cell):
        canonical = validate_and_canonicalize_operations(
            setting.symmetry_operations, cell
        )
        return build_context(
            canonical, cell,
            source_kind=SymmetrySourceKind.CATALOG_SETTING,
            setting_id=str(setting.setting_id),
            provenance=(("source", "space_group_setting"),),
        )
```

Normalize every operation to exact `(R, t mod 1)`, reject duplicates, validate integer rotations, determinant, identity, exact closure and inverses, and check `Rᵀ G R ≈ G` using one named metric tolerance. Sort by the exact descriptor before assigning keys and fingerprints. Do not add a public `SymmetryContext.status`; a successfully created context is valid, while provenance is represented by `source_kind`.

- [ ] **Step 4: Verify context tests and unchanged diffraction tests**

Run: `python3.11 -m pytest -q tests/orbit_first/test_symmetry_context.py`

Run: `python3.11 -m compileall -q src/cristma`

Expected: context tests pass and the complete package compiles. Run any diffraction tests present in the execution checkout explicitly; never mask a failing or absent path with `|| true`.

- [ ] **Step 5: Commit the context foundation**

```bash
git add src/cristma/crystallography/symmetry_context.py src/cristma/crystallography/__init__.py tests/orbit_first/test_symmetry_context.py
git commit -m "Add validated direct-space symmetry contexts"
```

### Task 2: Exact affine-periodic relation algebra

**Files:**
- Create: `src/cristma/crystallography/periodic_relation.py`
- Modify: `src/cristma/crystallography/__init__.py`
- Test: `tests/orbit_first/test_periodic_relation.py`

**Interfaces:**
- Consumes: `SymmetryContext.operation_keys` and its exact normalized operation lookup.
- Produces: `PeriodicSymmetryRelation`, `identity_relation(context)`, exact `compose`, `inverse`, `normalize`, and `apply` operations.

- [ ] **Step 1: Write failing group-law and lattice-carry tests**

```python
def test_relation_composition_rotates_lattice_translation_and_keeps_carry():
    left = PeriodicSymmetryRelation(rotation_key, (1, 0, 0))
    right = PeriodicSymmetryRelation(glide_key, (0, 1, 0))
    composed = left.compose(right, context)
    assert composed == hand_derived_relation

def test_relation_inverse_is_two_sided_exact_identity():
    inverse = relation.inverse(context)
    identity = identity_relation(context)
    assert relation.compose(inverse, context) == identity
    assert inverse.compose(relation, context) == identity
```

Use literal `Fraction` operations and a hand-derived nontrivial rotation plus fractional translation; do not calculate the expected relation with production helpers.

- [ ] **Step 2: Run the tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_periodic_relation.py`

Expected: import failure for `PeriodicSymmetryRelation`.

- [ ] **Step 3: Implement relation normalization, composition, and inversion**

```python
@dataclass(frozen=True, slots=True, order=True)
class PeriodicSymmetryRelation:
    operation_key: str
    lattice_translation: tuple[int, int, int]

    def compose(self, other, context) -> "PeriodicSymmetryRelation":
        return compose_periodic_relations(self, other, context)

    def inverse(self, context) -> "PeriodicSymmetryRelation":
        return invert_periodic_relation(self, context)

    def apply_fractional(self, coordinates, context) -> tuple[Fraction, Fraction, Fraction]:
        operation = context.operation_by_key(self.operation_key)
        return apply_operation_plus_lattice(
            operation, self.lattice_translation, coordinates
        )
```

Resolve both operation keys through the context. Compose exact affine actions,
normalize the fractional translation back to the context operation, and add
the exact integer carry. Reject relations referencing another context.

- [ ] **Step 4: Run exhaustive finite operation-pair checks**

Run: `python3.11 -m pytest -q tests/orbit_first/test_periodic_relation.py`

Add a parametrized check over every operation pair in representative catalog settings and verify closure, inverse, associativity samples, and integer translation output.

- [ ] **Step 5: Commit relation algebra**

```bash
git add src/cristma/crystallography/periodic_relation.py src/cristma/crystallography/__init__.py tests/orbit_first/test_periodic_relation.py
git commit -m "Add exact periodic symmetry relations"
```

### Task 3: Asymmetric-unit site mapping and stabilizers

**Files:**
- Create: `src/cristma/crystallography/asu_mapping.py`
- Modify: `src/cristma/crystallography/__init__.py`
- Test: `tests/orbit_first/test_asu_mapping.py`

**Interfaces:**
- Consumes: `CrystalStructure.sites`, `SymmetryContext`, and `PeriodicSymmetryRelation`.
- Produces: `SiteImage`, `SiteOrbitMapping`, `AsymmetricUnitMapping`, and `AsymmetricUnitMapper.build(structure, context)`.

- [ ] **Step 1: Write failing general- and special-position tests**

```python
def test_special_position_merges_images_but_retains_coset_evidence():
    mapping = AsymmetricUnitMapper().build(special_structure, context)
    site = mapping.by_site_id[special_site.id]
    assert len(site.reference_cell_images) == hand_expected_multiplicity
    assert len(site.stabilizer_relations) == hand_expected_stabilizer_order
    assert any(len(image.equivalent_relations) > 1 for image in site.reference_cell_images)

def test_mapping_identity_is_independent_of_operation_order():
    assert build(first_context).fingerprint == build(reordered_context).fingerprint
```

Include a relation that fixes a site modulo a nonzero lattice translation and assert the complete relation, not only its operation key, is retained in `stabilizer_relations`.

- [ ] **Step 2: Run the tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_asu_mapping.py`

Expected: imports fail because mapping types do not exist.

- [ ] **Step 3: Implement deterministic site-orbit mappings**

```python
@dataclass(frozen=True, slots=True)
class SiteImage:
    image_id: str
    representative_relation: PeriodicSymmetryRelation
    equivalent_relations: tuple[PeriodicSymmetryRelation, ...]
    fractional_position: tuple[float, float, float]
    normalization_translation: tuple[int, int, int]

@dataclass(frozen=True, slots=True)
class SiteOrbitMapping:
    independent_site_id: str
    stabilizer_relations: tuple[PeriodicSymmetryRelation, ...]
    reference_cell_images: tuple[SiteImage, ...]

@dataclass(frozen=True, slots=True)
class AsymmetricUnitMapping:
    site_orbits: tuple[SiteOrbitMapping, ...]
    symmetry_context_fingerprint: str
    fractional_tolerance: float
    fingerprint: str
```

Apply every exact relation representative, wrap numerically only at the final coordinate boundary, group periodically equal images with the explicit tolerance, and retain all exact relation evidence. Validate orbit-stabilizer counts.

- [ ] **Step 4: Verify all mapping invariants**

Run: `python3.11 -m pytest -q tests/orbit_first/test_asu_mapping.py`

Run a catalog sweep fixture that checks representative general and special Wyckoff positions for all 530 settings without requiring the pair finder.

- [ ] **Step 5: Commit ASU mapping**

```bash
git add src/cristma/crystallography/asu_mapping.py src/cristma/crystallography/__init__.py tests/orbit_first/test_asu_mapping.py
git commit -m "Map asymmetric-unit sites with exact stabilizers"
```

### Task 4: Complete cutoff-bounded ASU pair search

**Files:**
- Create: `src/cristma/crystallography/symmetry_pairs.py`
- Modify: `src/cristma/crystallography/__init__.py`
- Test: `tests/orbit_first/test_symmetry_pair_search.py`

**Interfaces:**
- Consumes: `CrystalStructure`, `SymmetryContext`, `AsymmetricUnitMapping`, cutoff, distance tolerance, and optional `max_candidates`.
- Produces: raw `SymmetryPairCandidate` records and complete/incomplete search provenance used by Task 5.

- [ ] **Step 1: Write failing skew-cell completeness and boundary tests**

```python
def test_skew_cell_buffer_matches_large_bruteforce_reference():
    found = SymmetryPairFinder(cutoff=3.1).find_candidates(structure, context, mapping)
    assert canonical_candidate_geometry(found.candidates) == hand_checked_bruteforce_pairs

def test_cutoff_boundary_uses_only_metric_tolerance():
    found = SymmetryPairFinder(cutoff=2.0, distance_tolerance=1e-10).find_candidates(
        boundary_structure, context, mapping
    )
    assert tuple(item.distance for item in found.candidates) == (2.0,)
```

The brute-force expected set must be created independently over a deliberately oversized translation range for a tiny triclinic fixture. Add a fixture requiring a translation outside `-1..1` to catch fixed-supercell implementations.

- [ ] **Step 2: Run the tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_symmetry_pair_search.py`

Expected: `SymmetryPairFinder` import fails.

- [ ] **Step 3: Implement mathematically complete translation bounds and bins**

```python
@dataclass(frozen=True, slots=True)
class SymmetryPairSearchPolicy:
    cutoff: float
    distance_tolerance: float = 1e-12
    max_candidates: int | None = None

@dataclass(frozen=True, slots=True)
class SymmetryPairCandidate:
    first_site_id: str
    second_site_id: str
    relation: PeriodicSymmetryRelation
    distance: float
    vector_cartesian: tuple[float, float, float]

@dataclass(frozen=True, slots=True)
class PairCandidateResult:
    candidates: tuple[SymmetryPairCandidate, ...]
    complete: bool
    integer_points_tested: int
    buffered_images: int
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]

class SymmetryPairFinder:
    def find_candidates(self, structure, context, mapping) -> PairCandidateResult:
        buffer = build_complete_metric_buffer(
            structure.cell, mapping, self.policy.cutoff,
            self.policy.distance_tolerance,
        )
        candidates = distance_checked_bin_pairs(buffer, self.policy)
        return PairCandidateResult(
            candidates=candidates.rows,
            complete=not candidates.limit_reached,
            integer_points_tested=candidates.integer_points_tested,
            buffered_images=len(buffer.images),
            diagnostics=candidates.diagnostics,
            provenance=candidates.provenance,
        )
```

Derive lattice-translation bounds from the inverse cell matrix and cutoff,
generate the complete finite image buffer, place Cartesian images into bins,
and run a final numerical metric check. `max_candidates` counts full candidate
pairs reaching the distance check. Reaching it returns controlled incomplete
provenance; it never returns a falsely complete result.

- [ ] **Step 4: Verify search completeness and benchmark candidate reduction**

Run: `python3.11 -m pytest -q tests/orbit_first/test_symmetry_pair_search.py`

Run a diagnostic benchmark on La2Zr2O7 and record independent-site count, buffered image count, candidate checks, and accepted raw relations. Do not yet compare final contact-orbit IDs because canonicalization is Task 5.

- [ ] **Step 5: Commit pair candidate search**

```bash
git add src/cristma/crystallography/symmetry_pairs.py src/cristma/crystallography/__init__.py tests/orbit_first/test_symmetry_pair_search.py
git commit -m "Search complete asymmetric-unit contact buffers"
```

### Task 5: Stabilizer quotient, ownership, and `SymmetryPairTable`

**Files:**
- Create: `src/cristma/crystallography/pair_canonical.py`
- Modify: `src/cristma/crystallography/symmetry_pairs.py`
- Modify: `src/cristma/crystallography/__init__.py`
- Test: `tests/orbit_first/test_pair_canonicalization.py`
- Test: `tests/orbit_first/test_symmetry_pair_table.py`

**Interfaces:**
- Consumes: `SymmetryPairCandidate`, endpoint `SiteOrbitMapping` stabilizer relations, and `SymmetryContext`.
- Produces: `canonical_pair_relation(...)`, `canonical_instance_owner(...)`, `SymmetryContactOrbit`, `SymmetryPairTable`, and `SymmetryPairFinder.find(...) -> SymmetryPairTable`.

- [ ] **Step 1: Write failing double-coset and endpoint-exchange tests**

```python
def test_special_position_descriptions_collapse_under_both_stabilizers():
    canonical = {
        canonical_pair_relation(candidate, first_mapping, second_mapping, context)
        for candidate in equivalent_raw_candidates
    }
    assert canonical == {hand_expected_descriptor}

def test_undirected_relation_and_inverse_have_one_identity():
    forward = canonical_pair_relation(ab_candidate, map_a, map_b, context)
    reverse = canonical_pair_relation(ba_inverse_candidate, map_b, map_a, context)
    assert forward == reverse
```

Use a hand-derived special-position fixture where different raw operations and nonzero lattice shifts describe the same physical pair.

- [ ] **Step 2: Run canonicalization tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_pair_canonicalization.py`

Expected: canonicalization functions do not exist.

- [ ] **Step 3: Implement exact quotient and ownership**

```python
def canonical_pair_relation(candidate, first_mapping, second_mapping, context):
    descriptors = []
    for left in first_mapping.stabilizer_relations:
        for right in second_mapping.stabilizer_relations:
            relation = left.compose(candidate.relation, context).compose(right, context)
            descriptors.extend((
                directed_descriptor(first_id, second_id, relation),
                directed_descriptor(second_id, first_id, relation.inverse(context)),
            ))
    return min(descriptors)
```

Implement the exact left/right action required by the chosen relation convention and verify it against the hand-derived fixture; do not copy the illustrative composition order blindly if the implemented source-to-target convention reverses it. Define one canonical global-translation and endpoint-exchange owner rule and expose it for both multiplicity counting and materialization.

- [ ] **Step 4: Write failing pair-table identity and multiplicity tests**

```python
def test_geometry_id_ignores_metric_only_cell_change():
    first = finder.find(structure_at_temperature_1, context_1, mapping_1)
    second = finder.find(structure_at_temperature_2, context_2, mapping_2)
    assert orbit_ids(first) == orbit_ids(second)
    assert distances(first) != distances(second)

def test_pair_multiplicity_is_not_chain_coordination_number():
    table = finder.find(one_site_periodic_chain, context, mapping)
    assert table.contact_orbits[0].multiplicity_in_reference_cell == 1
```

- [ ] **Step 5: Implement immutable geometric orbit and table models**

```python
class PairTableStatus(StrEnum):
    COMPLETE = "complete"
    INCOMPLETE = "incomplete"

@dataclass(frozen=True, slots=True)
class SymmetryContactOrbit:
    geometry_orbit_id: str
    first_independent_site_id: str
    second_independent_site_id: str
    canonical_relation: PeriodicSymmetryRelation
    equivalent_relations: tuple[PeriodicSymmetryRelation, ...]
    endpoint_stabilizers: tuple[tuple[PeriodicSymmetryRelation, ...], ...]
    representative_distance: float
    representative_vector_cartesian: tuple[float, float, float]
    multiplicity_in_reference_cell: int
    status: PairTableStatus
    diagnostics: tuple[Diagnostic, ...] = ()
    provenance: tuple[tuple[str, object], ...] = ()

@dataclass(frozen=True, slots=True)
class SymmetryPairTable:
    contact_orbits: tuple[SymmetryContactOrbit, ...]
    symmetry_context_fingerprint: str
    asymmetric_unit_mapping_fingerprint: str
    cutoff: float
    distance_tolerance: float
    status: PairTableStatus
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]
```

Hash IDs from canonical exact descriptors and the symmetry-action fingerprint. Exclude distance, cutoff, full context/cell fingerprint, and enumeration order.

Complete `SymmetryPairFinder.find` by passing `find_candidates` output through canonicalization and orbit aggregation:

```python
def find(self, structure, context, mapping) -> SymmetryPairTable:
    candidates = self.find_candidates(structure, context, mapping)
    return build_symmetry_pair_table(candidates, context, mapping, self.policy)
```

- [ ] **Step 6: Verify pair-table invariants and performance shape**

Run: `python3.11 -m pytest -q tests/orbit_first/test_pair_canonicalization.py tests/orbit_first/test_symmetry_pair_table.py`

For La2Zr2O7 assert that the pair table contains independent geometric orbits rather than 4096 expanded contacts and that no loop scales as `expanded_contact_count × operation_count × site_multiplicity`.

- [ ] **Step 7: Commit geometric pair tables**

```bash
git add src/cristma/crystallography/pair_canonical.py src/cristma/crystallography/symmetry_pairs.py src/cristma/crystallography/__init__.py tests/orbit_first/test_pair_canonicalization.py tests/orbit_first/test_symmetry_pair_table.py
git commit -m "Canonicalize asymmetric-unit contact orbits"
```

### Task 6: Orbit-first chemical interpretations

**Files:**
- Create: `src/cristma/crystal_chemistry/orbit_contacts.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Test: `tests/orbit_first/test_orbit_contact_resolution.py`

**Interfaces:**
- Consumes: `SymmetryPairTable`, `CrystalStructure` independent-site components, `CompositionGrammar`, `ReferenceData`, and `ShellResolutionPolicy` search limits.
- Produces: `OrientationMode`, `EndpointRole`, `ContactInterpretation`, `ResolvedContactOrbit`, and `ContactOrbitResolver.resolve(...)`.

- [ ] **Step 1: Write failing multi-interpretation and missing-radius tests**

```python
def test_mixed_pair_preserves_every_supported_interpretation():
    resolved = resolver.resolve(pair_table, mixed_structure, grammar)
    orbit = resolved.contact_orbits[0]
    assert tuple(item.interaction_type for item in orbit.interpretations) == (
        GrammarOperation.CENTRE_LIGAND_SHELL,
        GrammarOperation.COVALENT_NETWORK,
    )

def test_missing_radius_marks_only_affected_interpretation_incomplete():
    resolved = resolver.resolve(pair_table, structure_with_unknown_component, grammar)
    assert resolved.contact_orbits[0].status is ResolutionStatus.INCOMPLETE
    assert "crystal_chemistry.contact.radius_missing" in diagnostic_codes(resolved)
```

Expected interpretation order is literal and follows the canonical grammar scope order defined in the test fixture.

- [ ] **Step 2: Run the tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_orbit_contact_resolution.py`

Expected: orbit-contact models and resolver do not exist.

- [ ] **Step 3: Refactor component-pair interpretation to orbit input**

```python
class OrientationMode(StrEnum):
    UNDIRECTED = "undirected"
    ENDPOINT_ROLES = "endpoint_roles"

class EndpointRole(StrEnum):
    CENTER = "center"
    LIGAND = "ligand"
    NETWORK = "network"

@dataclass(frozen=True, slots=True)
class ContactInterpretation:
    interpretation_id: str
    interaction_type: GrammarOperation
    interaction_layer: InteractionLayer
    grammar_priority: InteractionPriority
    orientation_mode: OrientationMode
    endpoint_roles: tuple[EndpointRole, ...]
    component_pair_interpretations: tuple[ComponentPairInterpretation, ...]
    normalized_distance_range: tuple[float, float]
    status: ResolutionStatus
    evidence: tuple[SecondaryEvidence, ...]

@dataclass(frozen=True, slots=True)
class ResolvedContactOrbit:
    resolved_contact_orbit_id: str
    geometry_orbit_id: str
    interpretations: tuple[ContactInterpretation, ...]
    status: ResolutionStatus
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]
```

Move reusable chemistry from private expanded-contact helpers into the new resolver. Interpret each independent geometry orbit once. Preserve all distinct grammar scopes and component alternatives; do not call `_coalesce_resolved_contacts`.

- [ ] **Step 4: Verify chemistry without materialization**

Run: `python3.11 -m pytest -q tests/orbit_first/test_orbit_contact_resolution.py`

Add an assertion that no `ResolvedContact`, `GeometricContact`, or expanded atom ID occurs in the orbit result.

- [ ] **Step 5: Commit orbit chemical resolution**

```bash
git add src/cristma/crystal_chemistry/orbit_contacts.py src/cristma/crystal_chemistry/__init__.py tests/orbit_first/test_orbit_contact_resolution.py
git commit -m "Resolve contact chemistry by pair orbit"
```

### Task 7: Oriented contact incidence orbits

**Files:**
- Create: `src/cristma/crystal_chemistry/incidence_orbits.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Test: `tests/orbit_first/test_contact_incidences.py`

**Interfaces:**
- Consumes: `ResolvedContactOrbit`, its geometry orbit, `ContactInterpretation.endpoint_roles`, endpoint stabilizers, and independent-site occupancies.
- Produces: `ContactIncidenceOrbit` and `ContactIncidenceBuilder.build(...)`.

- [ ] **Step 1: Write failing chain, special-position, and mixed-occupancy tests**

```python
def test_one_pair_orbit_gives_two_local_incidences_in_simple_chain():
    incidences = builder.build(chain_pair_table, chain_contact_orbits, chain_structure)
    assert len(incidences) == 1
    assert incidences[0].incidence_multiplicity_per_center == 2

def test_effective_occupancy_uses_only_participating_component():
    incidence = only(builder.build(mixed_of_pair_table, interpretations, structure))
    assert incidence.effective_neighbor_occupancy == pytest.approx(0.6)

def test_incidence_identity_does_not_include_shell_policy():
    assert build(policy_a)[0].incidence_orbit_id == build(policy_b)[0].incidence_orbit_id
```

- [ ] **Step 2: Run the tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_contact_incidences.py`

Expected: incidence types do not exist.

- [ ] **Step 3: Implement centre-stabilizer incidence quotient**

```python
@dataclass(frozen=True, slots=True)
class ContactIncidenceOrbit:
    incidence_orbit_id: str
    resolved_contact_orbit_id: str
    interpretation_id: str
    center_independent_site_id: str
    ligand_independent_site_id: str
    oriented_periodic_relation: PeriodicSymmetryRelation
    incidence_multiplicity_per_center: int
    effective_neighbor_occupancy: float
    status: ResolutionStatus
    evidence: tuple[SecondaryEvidence, ...]
```

Orient each admissible interpretation around its scientific centre role, quotient neighbour images by the centre stabilizer, and count distinct local incidences. For undirected interpretations create the mathematically distinct endpoint-centred incidence orbits. Sum only participating ligand component occupancies; do not multiply by centre occupancy.

- [ ] **Step 4: Verify incidence identities and multiplicities**

Run: `python3.11 -m pytest -q tests/orbit_first/test_contact_incidences.py`

Run the special-position fixtures with reordered operation sets and assert identical incidence IDs, relations, and multiplicities.

- [ ] **Step 5: Commit incidence orbits**

```bash
git add src/cristma/crystal_chemistry/incidence_orbits.py src/cristma/crystal_chemistry/__init__.py tests/orbit_first/test_contact_incidences.py
git commit -m "Build oriented contact incidence orbits"
```

### Task 8: Shell alternatives and aggregate contact result

**Files:**
- Create: `src/cristma/crystal_chemistry/shell_orbits.py`
- Create: `src/cristma/crystal_chemistry/contact_analysis.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Test: `tests/orbit_first/test_shell_orbits.py`
- Test: `tests/orbit_first/test_contact_analysis_status.py`

**Interfaces:**
- Consumes: `SymmetryPairTable`, `ResolvedContactOrbit`, `ContactIncidenceOrbit`, and `ShellResolutionPolicy`.
- Produces: `ShellRole`, `CoordinationShellAlternative`, `CoordinationShellOrbit`, `ContactAnalysisResult`, `ContactAnalyzer.analyze(structure, symmetry_context, grammar)`.

- [ ] **Step 1: Write failing resolved/ambiguous weighted-shell tests**

```python
def test_ambiguous_shell_keeps_alternatives_without_selecting_one():
    shell = resolver.resolve(ambiguous_incidences)
    assert shell.status is ResolutionStatus.AMBIGUOUS
    assert shell.selected_alternative is None
    assert len(shell.alternatives) == 2
    assert same_incidence_id in shell.alternatives[0].primary_incidences
    assert same_incidence_id in shell.alternatives[1].secondary_incidences

def test_shell_cn_and_occupancy_use_incidence_weights():
    alternative = resolved_shell.selected
    assert alternative.geometric_CN == 4
    assert alternative.mean_occupied_neighbors == pytest.approx(3.2)
```

- [ ] **Step 2: Run shell tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_shell_orbits.py`

Expected: shell-orbit types do not exist.

- [ ] **Step 3: Implement multiplicity-weighted shell resolution**

```python
class ShellRole(StrEnum):
    PRIMARY = "primary"
    SECONDARY = "secondary"

@dataclass(frozen=True, slots=True)
class CoordinationShellAlternative:
    alternative_id: str
    primary_incidences: tuple[str, ...]
    secondary_incidences: tuple[str, ...]
    geometric_CN: int
    mean_occupied_neighbors: float
    boundary_evidence: tuple[SecondaryEvidence, ...]
    status: ResolutionStatus

@dataclass(frozen=True, slots=True)
class CoordinationShellOrbit:
    shell_orbit_id: str
    center_independent_site_id: str
    selected_alternative: str | None
    alternatives: tuple[CoordinationShellAlternative, ...]
    status: ResolutionStatus
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]

    @property
    def selected(self) -> CoordinationShellAlternative | None:
        if self.selected_alternative is None:
            return None
        return next(
            item for item in self.alternatives
            if item.alternative_id == self.selected_alternative
        )
```

Group unique incidence-orbit distances while retaining their local multiplicity as statistical weight. Move `PRIMARY` and `SECONDARY` assignment entirely into alternatives. Preserve lower/upper mixed-component alternatives rather than collapsing them.

- [ ] **Step 4: Write the failing top-level status truth table**

```python
@pytest.mark.parametrize((pair_status, shell_statuses, expected), (
    (PairTableStatus.INCOMPLETE, (), ResolutionStatus.INCOMPLETE),
    (PairTableStatus.COMPLETE, (ResolutionStatus.INCOMPLETE,), ResolutionStatus.INCOMPLETE),
    (PairTableStatus.COMPLETE, (ResolutionStatus.AMBIGUOUS,), ResolutionStatus.AMBIGUOUS),
    (PairTableStatus.COMPLETE, (ResolutionStatus.RESOLVED,), ResolutionStatus.RESOLVED),
    (PairTableStatus.COMPLETE, (), ResolutionStatus.NOT_APPLICABLE),
))
def test_contact_analysis_status_order(pair_status, shell_statuses, expected):
    assert aggregate_contact_analysis_status(pair_status, shell_statuses) is expected
```

- [ ] **Step 5: Implement `ContactAnalysisResult` and orchestration**

```python
@dataclass(frozen=True, slots=True)
class ContactAnalysisResult:
    _structure: CrystalStructure = field(repr=False, compare=False)
    _asymmetric_unit_mapping: AsymmetricUnitMapping = field(repr=False, compare=False)
    pair_table: SymmetryPairTable
    contact_orbits: tuple[ResolvedContactOrbit, ...]
    contact_incidence_orbits: tuple[ContactIncidenceOrbit, ...]
    coordination_shell_orbits: tuple[CoordinationShellOrbit, ...]
    status: ResolutionStatus
    diagnostics: tuple[Diagnostic, ...]
    configuration: tuple[tuple[str, object], ...]
    provenance: tuple[tuple[str, object], ...]

class ContactAnalyzer:
    def analyze(self, structure, symmetry_context, grammar) -> ContactAnalysisResult:
        mapping = AsymmetricUnitMapper().build(structure, symmetry_context)
        pair_table = self.pair_finder.find(structure, symmetry_context, mapping)
        contact_orbits = self.contact_resolver.resolve(pair_table, structure, grammar)
        incidences = self.incidence_builder.build(
            pair_table, contact_orbits, structure, mapping
        )
        shells = self.shell_resolver.resolve(pair_table, contact_orbits, incidences)
        return build_contact_analysis_result(
            structure, mapping, pair_table, contact_orbits, incidences, shells,
            configuration=self.get_config(),
        )
```

Aggregate status in this exact priority: incomplete pair table; any applicable incomplete shell; any applicable ambiguous shell; all applicable shells resolved; no applicable shells means `NOT_APPLICABLE`. Existing resolved contact orbits remain present when coordination is not applicable.

The two private immutable input references exist only so the compatibility property and explicit materializers can reconstruct coordinates and atom references. Validate their fingerprints against the pair table in `__post_init__`; exclude them from equality and public serialization. They are not a second calculated scientific state.

- [ ] **Step 6: Verify the complete orbit-first contact pipeline**

Run: `python3.11 -m pytest -q tests/orbit_first/test_shell_orbits.py tests/orbit_first/test_contact_analysis_status.py tests/orbit_first/test_orbit_contact_resolution.py tests/orbit_first/test_contact_incidences.py`

Assert that the analyzer takes an explicit context and that no expanded contact tuple is stored in the result.

- [ ] **Step 7: Commit shell and result orchestration**

```bash
git add src/cristma/crystal_chemistry/shell_orbits.py src/cristma/crystal_chemistry/contact_analysis.py src/cristma/crystal_chemistry/__init__.py tests/orbit_first/test_shell_orbits.py tests/orbit_first/test_contact_analysis_status.py
git commit -m "Resolve coordination shells by incidence orbit"
```

### Task 9: Contact and shell materializers

**Files:**
- Create: `src/cristma/crystal_chemistry/materialize.py`
- Modify: `src/cristma/crystal_chemistry/contact_analysis.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Test: `tests/orbit_first/test_contact_materialization.py`

**Interfaces:**
- Consumes: `ContactAnalysisResult`, `AsymmetricUnitMapping`, `ReferenceCell` or `CellRange`, and optional orbit/interpretation/alternative filters.
- Produces: `ReferenceCell`, `CellRange`, `ShellMembership`, materialized `ResolvedContact`, `materialize_contacts`, `materialize_coordination_shells`, and the read-only `ContactAnalysisResult.contacts` property.

- [ ] **Step 1: Write failing boundary-crossing ownership tests**

```python
def test_reference_cell_keeps_boundary_crossing_owned_contact_once():
    contacts = result.materialize_contacts(ReferenceCell())
    assert len(contacts) == orbit.multiplicity_in_reference_cell
    assert any(contact.second_atom_ref.cell_translation == (1, 0, 0) for contact in contacts)

def test_larger_region_translates_owner_cells_without_duplicates():
    contacts = result.materialize_contacts(CellRange(0, 2, 0, 1, 0, 1))
    assert len({contact.contact_id for contact in contacts}) == len(contacts)
    assert len(contacts) == 2 * orbit.multiplicity_in_reference_cell
```

- [ ] **Step 2: Run materialization tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_contact_materialization.py`

Expected: materializer and region types do not exist.

- [ ] **Step 3: Implement outward-only deterministic materialization**

```python
@dataclass(frozen=True, slots=True)
class ShellMembership:
    shell_orbit_id: str
    alternative_id: str
    role: ShellRole

@dataclass(frozen=True, slots=True)
class ReferenceCell:
    pass

@dataclass(frozen=True, slots=True)
class CellRange:
    a_min: int
    a_stop: int
    b_min: int
    b_stop: int
    c_min: int
    c_stop: int

@dataclass(frozen=True, slots=True)
class ResolvedContact:
    contact_id: str
    resolved_contact_orbit_id: str
    first_atom_ref: PeriodicAtomRef
    second_atom_ref: PeriodicAtomRef
    distance: float
    vector_cartesian: tuple[float, float, float]
    interpretations: tuple[ContactInterpretation, ...]
    shell_memberships: tuple[ShellMembership, ...]
    provenance: tuple[tuple[str, object], ...]

def materialize_contacts(
    result,
    structure,
    mapping,
    region=ReferenceCell(),
    **filters,
):
    owned = owned_pair_instances(result.pair_table, mapping, region, filters)
    return tuple(
        materialize_resolved_instance(item, result, structure, mapping)
        for item in owned
    )
```

Validate `CellRange` as nonempty half-open integer intervals. Thus
`CellRange(0, 2, 0, 1, 0, 1)` contains exactly two owner cells. Implement
`materialize_coordination_shells` by grouping the already materialized contact
instances through their `ShellMembership` records; it must not recalculate a
shell boundary:

```python
def materialize_coordination_shells(result, region=ReferenceCell(), **filters):
    contacts = result.materialize_contacts(region, **filters)
    return group_materialized_shell_instances(
        contacts,
        result.coordination_shell_orbits,
        region,
    )
```

Use Task 5's ownership function; do not independently reimplement symmetry canonicalization. Permit one endpoint outside the owner region. Hash instance IDs from resolved orbit ID plus canonical periodic endpoint references and sort output deterministically.

- [ ] **Step 4: Add and verify the compatibility property**

```python
@property
def contacts(self) -> tuple[ResolvedContact, ...]:
    return materialize_contacts(
        self,
        self._structure,
        self._asymmetric_unit_mapping,
        ReferenceCell(),
    )

def materialize_contacts(self, region=ReferenceCell(), **filters):
    return materialize_contacts(
        self,
        self._structure,
        self._asymmetric_unit_mapping,
        region,
        **filters,
    )
```

Store the mapping or an immutable reference sufficient for materialization in the result orchestration without storing materialized contacts. Test `result.contacts == result.materialize_contacts(ReferenceCell())` and inspect dataclass fields to prove no `contacts` field exists.

- [ ] **Step 5: Verify materialization invariants**

Run: `python3.11 -m pytest -q tests/orbit_first/test_contact_materialization.py`

Check reference-cell count equals every orbit multiplicity, all IDs are unique, ambiguous alternatives produce distinct shell memberships, and materialization results never appear in analyzer inputs.

- [ ] **Step 6: Commit materializers**

```bash
git add src/cristma/crystal_chemistry/materialize.py src/cristma/crystal_chemistry/contact_analysis.py src/cristma/crystal_chemistry/__init__.py tests/orbit_first/test_contact_materialization.py
git commit -m "Materialize contact orbits for consumers"
```

### Task 10: Orbit-first coordination polyhedra

**Files:**
- Modify: `src/cristma/crystal_chemistry/polyhedra.py`
- Modify: `src/cristma/crystal_chemistry/polyhedron_orbits.py`
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Test: `tests/orbit_first/test_orbit_polyhedra.py`

**Interfaces:**
- Consumes: `CoordinationShellOrbit`, its selected alternative, `ContactIncidenceOrbit`, `AsymmetricUnitMapping`, and the immutable source structure.
- Produces: representative `CoordinationPolyhedron`, `CoordinationPolyhedronOrbit`, and `PolyhedronOrbitBuilder.build(contact_result)` without materialized global contacts.

- [ ] **Step 1: Write failing local-realization and ambiguity tests**

```python
def test_polyhedron_realizes_local_vertices_without_global_contacts():
    result = PolyhedronOrbitBuilder().build(resolved_contact_result)
    polyhedron = result.polyhedron_orbits[0].representative
    assert polyhedron.coordination_number == expected_cn
    assert tuple(vertex.atom_ref for vertex in polyhedron.vertices) == expected_local_refs
    assert not hasattr(polyhedron, "vertex_contacts")

def test_ambiguous_shell_does_not_guess_one_polyhedron():
    result = PolyhedronOrbitBuilder().build(ambiguous_contact_result)
    assert result.polyhedron_orbits == ()
    assert "crystal_chemistry.polyhedron.shell_ambiguous" in result.diagnostic_codes
```

- [ ] **Step 2: Run tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_orbit_polyhedra.py`

Expected: current builder requires expanded `CoordinationShell` and cannot consume shell orbits.

- [ ] **Step 3: Replace expanded-contact vertex construction with local incidence realization**

```python
class PolyhedronOrbitBuilder:
    def build(self, contact_result: ContactAnalysisResult) -> PolyhedronOrbitBuildResult:
        for shell in contact_result.coordination_shell_orbits:
            if shell.status is not ResolutionStatus.RESOLVED:
                continue
            alternative = shell.selected
            vertices = realize_local_primary_incidences(
                shell, alternative, contact_result.contact_incidence_orbits,
                contact_result._asymmetric_unit_mapping,
            )
            polyhedron = build_polyhedron_geometry(shell, alternative, vertices)
            representative_rows.append(polyhedron)
        return group_representatives_by_shell_orbit(representative_rows)
```

Realize only the finite local neighbours around one representative centre. Preserve existing convex-hull faces, Baur distortion, edge-angle dispersion, volume, centroid, centre offset, mixed/split diagnostics, and canonical face signature. Replace `vertex_contacts` with incidence-orbit references and local `PeriodicAtomRef` values.

- [ ] **Step 4: Verify descriptor parity where legacy is known correct**

Run: `python3.11 -m pytest -q tests/orbit_first/test_orbit_polyhedra.py`

Compare resolved non-special fixtures against hand-calculated vertices and existing descriptors. For special positions, validate analytical CN and unique local vertices rather than forcing legacy equality.

- [ ] **Step 5: Commit orbit-first polyhedra**

```bash
git add src/cristma/crystal_chemistry/polyhedra.py src/cristma/crystal_chemistry/polyhedron_orbits.py src/cristma/crystal_chemistry/__init__.py tests/orbit_first/test_orbit_polyhedra.py
git commit -m "Build coordination polyhedra from shell orbits"
```

### Task 11: Orbit-first structural units, graph, and connectivity

**Files:**
- Modify: `src/cristma/crystal_chemistry/structural_units.py`
- Modify: `src/cristma/crystal_chemistry/structural_graph.py`
- Modify: `src/cristma/crystal_chemistry/representation.py`
- Modify: `src/cristma/crystal_chemistry/periodic_connectivity.py`
- Test: `tests/orbit_first/test_orbit_structural_graph.py`

**Interfaces:**
- Consumes: `ContactAnalysisResult`, `CoordinationPolyhedronOrbit`, geometric pair relations, and selected interaction layers.
- Produces: structural-unit orbits, an orbit-referenced structural graph, selected representations, and periodic components/ranks without materialized contacts.

- [ ] **Step 1: Write failing unit and translation-rank tests**

```python
def test_structural_graph_edges_reference_contact_orbits_only():
    graph = StructuralGraphBuilder().build(contact_result, polyhedron_result)
    assert all(edge.resolved_contact_orbit_id for edge in graph.connections)
    assert all(not hasattr(edge, "contact_id") for edge in graph.connections)

@pytest.mark.parametrize((fixture, expected_rank), (
    (finite_cluster, 0),
    (periodic_chain, 1),
    (periodic_layer, 2),
    (periodic_framework, 3),
))
def test_periodic_rank_comes_from_exact_relation_translations(fixture, expected_rank):
    assert analyze(fixture).components[0].rank == expected_rank
```

- [ ] **Step 2: Run hierarchy-core tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_orbit_structural_graph.py`

Expected: current builders require expanded resolved contacts and units.

- [ ] **Step 3: Migrate structural units and connections to orbit references**

```python
@dataclass(frozen=True, slots=True)
class StructuralConnection:
    first_unit_orbit_id: str
    second_unit_orbit_id: str
    resolved_contact_orbit_id: str
    relation: PeriodicSymmetryRelation
    kind: StructuralConnectionKind
```

Build unit identities from shell/polyhedron/contact orbit identities and exact relations. Keep planar BO3 and polyhedral BO4 geometry scientifically determined from local coordinates, but do not create symmetry-expanded unit collections.

- [ ] **Step 4: Derive periodic connectivity from exact cycle translations**

```python
translation_generators = exact_cycle_translation_generators(component)
rank = integer_lattice_rank(translation_generators)
```

Use integer relation translations and graph cycles; numerical Cartesian coordinates must not determine topological rank. Preserve simultaneous child motifs inside a rank-3 parent for the later motif layer without enumerating arbitrary subgraphs.

- [ ] **Step 5: Verify unit orbits and ranks**

Run: `python3.11 -m pytest -q tests/orbit_first/test_orbit_structural_graph.py`

Check BO3/BO4 unit counts and rank classifications on LiB3O5, K7, natrolite, NaCl, and synthetic finite/chain/layer/framework fixtures.

- [ ] **Step 6: Commit orbit-first hierarchy core**

```bash
git add src/cristma/crystal_chemistry/structural_units.py src/cristma/crystal_chemistry/structural_graph.py src/cristma/crystal_chemistry/representation.py src/cristma/crystal_chemistry/periodic_connectivity.py tests/orbit_first/test_orbit_structural_graph.py
git commit -m "Build structural graphs from scientific orbits"
```

### Task 12: Orbit-first blocks and rings

**Files:**
- Modify: `src/cristma/crystal_chemistry/structural_blocks.py`
- Modify: `src/cristma/crystal_chemistry/ring_finder.py`
- Modify: `src/cristma/crystal_chemistry/_ring_search.py`
- Modify: `src/cristma/crystal_chemistry/_ring_symmetry.py`
- Modify: `src/cristma/crystal_chemistry/rings.py`
- Test: `tests/orbit_first/test_orbit_blocks_and_rings.py`

**Interfaces:**
- Consumes: orbit-first structural graph, selected representation, periodic components, and exact relations.
- Produces: block orbits and ring orbits whose parent and member identities are already symmetry-native.

- [ ] **Step 1: Write failing scientific ring regressions**

```python
def test_nacl_coordination_packing_has_no_structural_rings():
    assert analyze_structure(nacl).ring_orbits == ()

def test_k7_keeps_one_isolated_borate_ring_orbit():
    result = analyze_structure(k7)
    assert len(result.rings) == 18
    assert len(result.ring_orbits) == 1

def test_natrolite_framework_ring_analysis_is_preserved():
    result = analyze_structure(natrolite)
    assert result.periodic_rank == 3
    assert result.ring_orbits
```

- [ ] **Step 2: Run block/ring tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_orbit_blocks_and_rings.py`

Expected: existing finders require expanded blocks, units, or contacts.

- [ ] **Step 3: Migrate block and ring identity to orbit relations**

```python
blocks = StructuralBlockFinder().find(representation, connectivity)
rings = RingFinder().find(representation, blocks)
```

Remove atomic-view and expanded-structure arguments from scientific ring identity. Search cycles in the structural orbit graph, retaining exact accumulated relation translations. Keep pure coordination edges excluded from structural rings. Build `parent_block_orbit_id` directly rather than grouping expanded rings after search.

- [ ] **Step 4: Verify deterministic rings and no combinatorial regression**

Run: `python3.11 -m pytest -q tests/orbit_first/test_orbit_blocks_and_rings.py`

Assert ring IDs and ordering survive operation reordering. Retain the acyclic precheck. Benchmark NaCl, K7, and natrolite and ensure orbit-first traversal does not materialize all symmetry-expanded rings during scientific analysis.

- [ ] **Step 5: Commit orbit-first blocks and rings**

```bash
git add src/cristma/crystal_chemistry/structural_blocks.py src/cristma/crystal_chemistry/ring_finder.py src/cristma/crystal_chemistry/_ring_search.py src/cristma/crystal_chemistry/_ring_symmetry.py src/cristma/crystal_chemistry/rings.py tests/orbit_first/test_orbit_blocks_and_rings.py
git commit -m "Analyze structural blocks and rings by orbit"
```

### Task 13: Public cutover and legacy deletion

**Files:**
- Modify: `src/cristma/crystal_chemistry/__init__.py`
- Modify: `src/cristma/crystal_chemistry/contacts.py` or delete it after moving retained value types.
- Delete: `src/cristma/crystal_chemistry/resolver.py`
- Delete: `src/cristma/crystal_chemistry/contact_orbits.py`
- Modify: `tools/smoke_cif_corpus.py`
- Modify: `README.md`
- Test: `tests/orbit_first/test_public_orbit_api.py`
- Test: `tests/orbit_first/test_no_expanded_scientific_route.py`

**Interfaces:**
- Consumes: all orbit-first components from Tasks 1–12.
- Produces: the final public `ContactAnalyzer -> ContactAnalysisResult` route and compatibility materializers; removes the second scientific API.

- [ ] **Step 1: Write failing public API and architecture tests**

```python
def test_public_contact_route_is_orbit_first():
    result = ContactAnalyzer(policy).analyze(structure, symmetry_context, grammar)
    assert result.contact_orbits
    assert result.contacts == result.materialize_contacts(ReferenceCell())

def test_legacy_expanded_first_symbols_are_not_public():
    assert not hasattr(cristma.crystal_chemistry, "CoordinationShellResolver")
    assert not hasattr(cristma.crystal_chemistry, "CrystalChemistryResolution")
```

The architecture test must import and exercise the production API. Do not merely grep source text. A separate static import-boundary check may supplement it to ensure scientific modules do not import materialized `ResolvedContact`.

- [ ] **Step 2: Run cutover tests and verify RED**

Run: `python3.11 -m pytest -q tests/orbit_first/test_public_orbit_api.py tests/orbit_first/test_no_expanded_scientific_route.py`

Expected: legacy symbols still exist and production smoke route still uses them.

- [ ] **Step 3: Switch exports and smoke route in one change**

```python
context = SymmetryContext.from_definition(structure.space_group, structure.cell)
result = ContactAnalyzer(policy).analyze(structure, context, chemistry.grammar)
polyhedra = PolyhedronOrbitBuilder().build(result)
```

Update public exports and the corpus smoke tool. Update README examples to show explicit context creation, orbit-first analysis, and outward-only materialization. Do not add consumer-specific migration code.

- [ ] **Step 4: Delete legacy expanded-first code**

Remove `resolver.py`, `contact_orbits.py`, `_coalesce_resolved_contacts`, expanded shell construction, and old stored-contact result models. Move only still-valid small value types such as `ComponentPairInterpretation`, `SecondaryEvidence`, and `ResolutionStatus` into their orbit-first owner modules. Remove all scientific imports of materialized `ResolvedContact`.

- [ ] **Step 5: Run public and architecture tests**

Run: `python3.11 -m pytest -q tests/orbit_first/test_public_orbit_api.py tests/orbit_first/test_no_expanded_scientific_route.py`

Run: `python3.11 -m compileall -q src/cristma tools`

Expected: public route passes, legacy imports fail as asserted, and compilation succeeds.

- [ ] **Step 6: Commit the atomic cutover**

```bash
git add -A src/cristma/crystal_chemistry tools/smoke_cif_corpus.py README.md tests/orbit_first
git commit -m "Cut crystal chemistry over to orbit-first analysis"
```

### Task 14: Scientific corpus and performance acceptance

**Files:**
- Create: `tools/benchmark_orbit_contacts.py`
- Test: `tests/orbit_first/test_real_structure_acceptance.py`
- Modify: `README.md` only if measured behavior requires correcting an already written example.

**Interfaces:**
- Consumes: final public API only.
- Produces: reproducible acceptance evidence and a benchmark report; no new scientific model.

- [ ] **Step 1: Write real-structure acceptance assertions**

```python
@pytest.mark.parametrize("path", required_real_structure_paths())
def test_required_structure_completes_or_reports_scientific_status(path):
    result = analyze_path(path)
    assert result.status in set(ResolutionStatus)
    assert all(orbit.geometry_orbit_id for orbit in result.contact_orbits)

def test_real_ring_expectations():
    assert ring_summary(k7) == {"rings": 18, "orbits": 1}
    assert ring_summary(nacl) == {"rings": 0, "orbits": 0}
    assert ring_summary(natrolite)["orbits"] > 0
```

Include La2Zr2O7, La0.5Zr0.5O1.75, LiB3O5, Na6Mo11O36, NaLiMoO, grossular, gehlenite, K7, NaCl, natrolite, organic, mixed/split occupancy, and special-position fixtures. Missing optional external corpus paths must be explicitly skipped with their path in the reason; bundled analytical fixtures must never skip.

- [ ] **Step 2: Run analytical and curated acceptance tests**

Run: `python3.11 -m pytest -q tests/orbit_first`

Expected: every orbit-first invariant and available curated structure passes. Diagnose differences from legacy scientifically; update expectations only with a written analytical reason in the test.

- [ ] **Step 3: Implement the reproducible benchmark tool**

```python
for path in arguments.paths:
    started = time.perf_counter()
    structure, result = analyze_path(path)
    elapsed = time.perf_counter() - started
    print(json.dumps({
        "file": path.name,
        "independent_sites": len(structure.sites),
        "pair_orbits": len(result.pair_table.contact_orbits),
        "contact_orbits": len(result.contact_orbits),
        "seconds": elapsed,
    }, sort_keys=True))
```

Keep timing outside scientific provenance. Report Python version, platform, CrIStMa commit, policy, and input path/checksum. Do not enforce one absolute CI time across machines.

- [ ] **Step 4: Run La2Zr2O7 and hierarchy benchmarks**

Run: `PYTHONPATH=src python3.11 tools/benchmark_orbit_contacts.py '/Users/artem/Desktop/Cif/-_La2_O7_Zr2_-.cif' '/Users/artem/Desktop/Cif/K7CaY2(B5O10)3.cif' '/Users/artem/Desktop/Cif/Natrolite_deuterated_.cif'`

Expected: La2Zr2O7 contact analysis is at least one order of magnitude faster than the measured 85-second instrumented baseline on the same machine, while K7 and natrolite scientific summaries remain correct.

- [ ] **Step 5: Run both maintained external corpora**

Run: `PYTHONPATH=src python3.11 tools/smoke_cif_corpus.py '/Users/artem/Desktop/Cif'`

Run: `PYTHONPATH=src python3.11 tools/smoke_cif_corpus.py '/Users/artem/Library/Application Support/Sci/apps/xrd_phase_finder/data/cod_cache/cif' --skip-rings`

Expected: all readable structures remain available, intentionally empty coordinate sets remain nonfatal, and each failure is reported by exact file and diagnostic rather than stopping the corpus.

- [ ] **Step 6: Run final package verification**

Run: `python3.11 -m compileall -q src/cristma tools`

Run: `git diff --check`

Build an isolated wheel and sdist using the established beta build process, inspect their file lists, install the wheel into a clean environment, and execute the public orbit-first smoke example. Do not upload or change the package version in this task.

- [ ] **Step 7: Commit benchmark and acceptance coverage**

```bash
git add tools/benchmark_orbit_contacts.py tests/orbit_first/test_real_structure_acceptance.py README.md
git commit -m "Verify orbit-first crystal chemistry acceptance"
```

## Final Completion Gate

Before calling the migration complete, verify all of the following in one fresh run:

- every test under `tests/orbit_first` passes;
- catalog setting order invariance and special-position stabilizer tests pass;
- exact relation group-law tests pass;
- pair ownership and materialization-count invariants pass;
- no scientific module consumes materialized contacts;
- old public resolver/result/orbit-grouping symbols are absent;
- K7 has its accepted isolated ring result and NaCl has no structural rings;
- La2Zr2O7 meets the order-of-magnitude performance target;
- both available CIF corpora complete with only explicit per-file diagnostics;
- wheel and sdist contain runtime files only and install in a clean environment;
- the working tree is clean and the implementation commits are reviewable in task order.
