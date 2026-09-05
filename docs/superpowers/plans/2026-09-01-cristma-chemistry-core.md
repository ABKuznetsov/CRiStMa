# CrIStMa Chemistry Core Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a self-contained Chemistry layer that turns canonical structure composition into one actionable material family and concrete candidate interactions without inspecting coordinates.

**Architecture:** Immutable reference catalogs supply element categories, covalent radii, and the validated CRAFT v3/v3.1 knowledge base. Stateless Chemistry functions build a normalized `Composition`, classify it, and compile a composition grammar; `ChemistryAnalyzer` only composes those functions.

**Tech Stack:** Python 3.11+, dataclasses, `importlib.resources`, JSON, pytest, setuptools; NumPy remains the only runtime dependency.

**Spec:** `docs/superpowers/specs/2026-09-01-cristma-chemistry-and-reference-data-design.md`

## Global Constraints

- Chemistry never reads coordinates, distances, neighbors, bonds, polyhedra, or topology.
- Reference Data is immutable inert data and performs no runtime network access.
- No Qt, pymatgen, CRAFT, or application imports.
- Ordinary compositions return one actionable primary family; alternatives exist only for explicit curated grey-zone rules.
- Missing reference values fail explicitly; no guessed radius fallback.
- Display colors remain outside CrIStMa.

---

### Task 1: Element and covalent-radius reference catalogs

**Files:**
- Create: `src/cristma/reference_data/__init__.py`
- Create: `src/cristma/reference_data/elements.py`
- Create: `src/cristma/reference_data/radii.py`
- Create: `src/cristma/reference_data/resources/elements.json`
- Create: `src/cristma/reference_data/resources/covalent_radii.json`
- Test: `tests/reference_data/test_elements.py`
- Test: `tests/reference_data/test_radii.py`

**Interfaces:**
- Produces: `ElementCategory`, `ElementRecord`, `ElementCatalog.by_symbol(symbol)`, `CovalentRadiusRecord`, `CovalentRadii.find(symbol)`.
- `ElementRecord.is_metal` is the only metal/nonmetal policy consumed by classification and grammar.

- [ ] **Step 1: Write failing catalog tests**

```python
def test_element_catalog_identifies_fe_and_si():
    catalog = ElementCatalog.default()
    assert catalog.by_symbol("Fe").is_metal
    assert not catalog.by_symbol("Si").is_metal
    assert catalog.by_symbol("fe").symbol == "Fe"

def test_covalent_radius_has_units_and_no_guess():
    radii = CovalentRadii.default()
    assert radii.find("O").value == 0.66
    assert radii.find("O").unit == "angstrom"
    with pytest.raises(KeyError):
        radii.find("Xx")
```

- [ ] **Step 2: Run the focused tests and verify missing-module failures**

Run: `pytest tests/reference_data/test_elements.py tests/reference_data/test_radii.py -q`

- [ ] **Step 3: Implement frozen records and resource-backed catalogs**

```python
class ElementCategory(StrEnum):
    METAL = "metal"
    METALLOID = "metalloid"
    NONMETAL = "nonmetal"
    NOBLE_GAS = "noble_gas"

@dataclass(frozen=True, slots=True)
class ElementRecord:
    symbol: str
    atomic_number: int
    category: ElementCategory
    dataset_id: str
    dataset_version: str

    @property
    def is_metal(self) -> bool:
        return self.category is ElementCategory.METAL
```

- [ ] **Step 4: Run focused tests**

Run: `pytest tests/reference_data/test_elements.py tests/reference_data/test_radii.py -q`

- [ ] **Step 5: Commit**

```bash
git add src/cristma/reference_data tests/reference_data
git commit -m "feat(reference): add element and radius catalogs"
```

---

### Task 2: Validated Chemical Reference DB

**Files:**
- Create: `src/cristma/reference_data/chemical_reference.py`
- Create: `src/cristma/reference_data/facade.py`
- Create: `src/cristma/reference_data/resources/chemical_reference_v3.json`
- Create: `src/cristma/reference_data/resources/chemical_reference_v3_1.json`
- Modify: `src/cristma/reference_data/__init__.py`
- Test: `tests/reference_data/test_chemical_reference.py`

**Interfaces:**
- Produces: `ChemicalReference`, `ChemicalReferenceIntegrityReport`, `load_chemical_reference(version)`, `validate_reference_integrity(reference)`, `ReferenceData.default()`.
- `ReferenceData` exposes `elements`, `covalent_radii`, and `chemical` only.

- [ ] **Step 1: Write failing integrity and facade tests**

```python
def test_default_reference_is_validated_v31():
    reference = ReferenceData.default()
    assert reference.chemical.schema_version == "3.1.0-draft"
    assert reference.chemical.family("inorganic.oxide")["profile_id"]
    assert reference.chemical.boundary_case("CaSi2")["refined"]["preferred_candidates"] == (
        "inorganic.tetrelide", "inorganic.zintl"
    )

def test_historical_v3_counts_are_stable():
    report = validate_reference_integrity(load_chemical_reference("3.0.0-draft"))
    assert (report.family_count, report.group_count, report.boundary_case_count) == (103, 243, 155)
```

- [ ] **Step 2: Run the focused test and verify failure**

Run: `pytest tests/reference_data/test_chemical_reference.py -q`

- [ ] **Step 3: Migrate the two inert resources and loader, correcting only `CaSi2` in the supplied draft**

The loader freezes nested mappings/lists, validates schema `3.x`, family and group parent acyclicity, profile routes, boundary references, grammar selectors, operations, routes, and templates. It computes SHA-256 from packaged bytes.

- [ ] **Step 4: Implement the immutable facade**

```python
@dataclass(frozen=True, slots=True)
class ReferenceData:
    elements: ElementCatalog
    covalent_radii: CovalentRadii
    chemical: ChemicalReference

    @classmethod
    @lru_cache(maxsize=1)
    def default(cls) -> "ReferenceData":
        return cls(ElementCatalog.default(), CovalentRadii.default(), load_chemical_reference())
```

- [ ] **Step 5: Run focused reference tests and commit**

Run: `pytest tests/reference_data -q`

```bash
git add src/cristma/reference_data tests/reference_data
git commit -m "feat(reference): add chemical reference database"
```

---

### Task 3: Immutable composition and structure adapter

**Files:**
- Create: `src/cristma/chemistry/composition.py`
- Modify: `src/cristma/chemistry/__init__.py`
- Test: `tests/chemistry/test_composition.py`

**Interfaces:**
- Consumes: existing `ChemicalSpecies`/`SiteComponent` structural contracts and `ElementCatalog`.
- Produces: `Composition.from_mapping(values)` and `Composition.from_structure(structure)` with `amount(symbol)`, `elements`, and `normalized_formula`.

- [ ] **Step 1: Write failing normalization and multiplicity tests**

```python
def test_composition_normalizes_symbols_and_formula():
    composition = Composition.from_mapping({"ca": 1, "O": 1})
    assert composition.elements == ("Ca", "O")
    assert composition.normalized_formula == "CaO"

def test_independent_sites_use_occupancy_times_calculated_multiplicity(crystal):
    composition = Composition.from_structure(crystal)
    assert composition.amount("Ca") == pytest.approx(2.0)

def test_unknown_occupied_species_is_rejected(structure_with_unknown):
    with pytest.raises(ValueError, match="no known element"):
        Composition.from_structure(structure_with_unknown)
```

- [ ] **Step 2: Run the focused tests and verify failure**

Run: `pytest tests/chemistry/test_composition.py -q`

- [ ] **Step 3: Implement canonical immutable amounts and structural typing adapter**

The adapter reads only `components`, component `occupancy`, species element identity, and site multiplicity. Molecular atoms use multiplicity one. Independent sites prefer `calculated_multiplicity`; identity-only crystals use one; nonidentity sites without calculated multiplicity raise `ValueError`.

- [ ] **Step 4: Run tests and commit**

Run: `pytest tests/chemistry/test_composition.py -q`

```bash
git add src/cristma/chemistry/composition.py src/cristma/chemistry/__init__.py tests/chemistry/test_composition.py
git commit -m "feat(chemistry): add canonical composition"
```

---

### Task 4: Actionable material classification

**Files:**
- Create: `src/cristma/chemistry/evidence.py`
- Create: `src/cristma/chemistry/classification.py`
- Modify: `src/cristma/chemistry/__init__.py`
- Test: `tests/chemistry/test_classification.py`

**Interfaces:**
- Consumes: `Composition`, `ReferenceData`.
- Produces: `ChemicalEvidence`, `CompositionKind`, `ChemicalDomain`, `ChemicalClassification`, `classify_composition(composition, reference)`.
- `ChemicalClassification.primary_family` is a validated string ID; `alternative_families` is empty except for an explicit boundary rule.

- [ ] **Step 1: Write exact failing acceptance tests**

```python
@pytest.mark.parametrize(("formula", "family"), [
    ({"Fe": 1}, "elemental.metallic"),
    ({"Si": 1}, "elemental.covalent"),
    ({"Fe": 1, "Al": 1}, "inorganic.intermetallic"),
    ({"Ca": 1, "O": 1}, "inorganic.oxide"),
    ({"Na": 1, "Cl": 1}, "inorganic.halide"),
    ({"Fe": 1, "S": 2}, "inorganic.chalcogenide"),
])
def test_primary_material_family_is_actionable(formula, family):
    result = classify_composition(Composition.from_mapping(formula), ReferenceData.default())
    assert result.primary_family == family
    assert result.alternative_families == ()
```

- [ ] **Step 2: Run the focused test and verify failure**

Run: `pytest tests/chemistry/test_classification.py -q`

- [ ] **Step 3: Implement deterministic ordered rules backed by reference element sets**

Rule order is elemental identity, explicit curated grey-zone match, oxygen, halogen, nitrogen, other chalcogen, other pnictogen, carbon inorganic/organic signatures, all-metal intermetallic, then unresolved. Every result has evidence and method/reference provenance; limitations use the existing `Diagnostic` type.

- [ ] **Step 4: Run tests and commit**

Run: `pytest tests/chemistry/test_classification.py -q`

```bash
git add src/cristma/chemistry tests/chemistry/test_classification.py
git commit -m "feat(chemistry): classify material families"
```

---

### Task 5: Concrete interaction grammar and analyzer

**Files:**
- Create: `src/cristma/chemistry/grammar.py`
- Create: `src/cristma/chemistry/analyzer.py`
- Modify: `src/cristma/chemistry/__init__.py`
- Modify: `src/cristma/__init__.py`
- Test: `tests/chemistry/test_grammar.py`
- Test: `tests/chemistry/test_analyzer.py`

**Interfaces:**
- Consumes: `Composition`, `ChemicalClassification`, `ReferenceData`.
- Produces: `CandidateInteraction`, `CompositionGrammar`, `compile_composition_grammar`, `ChemistryResult`, `ChemistryAnalyzer.analyze(composition)`.

- [ ] **Step 1: Write exact interaction tests**

```python
@pytest.mark.parametrize(("formula", "pairs"), [
    ({"Fe": 1}, {(("Fe",), ("Fe",), "metallic_coordination", "primary")}),
    ({"Si": 1}, {(("Si",), ("Si",), "covalent_network", "primary")}),
    ({"Fe": 1, "Al": 1}, {(("Al",), ("Fe",), "metallic_coordination", "primary")}),
    ({"Ca": 1, "O": 1}, {(("Ca",), ("O",), "centre_ligand_shell", "primary")}),
    ({"Na": 1, "Cl": 1}, {(("Cl",), ("Na",), "centre_ligand_shell", "primary")}),
    ({"Fe": 1, "S": 2}, {
        (("Fe",), ("S",), "centre_ligand_shell", "primary"),
        (("S",), ("S",), "intra_subsystem_bonds", "allowed"),
    }),
])
def test_grammar_returns_concrete_search_pairs(formula, pairs):
    result = ChemistryAnalyzer().analyze(Composition.from_mapping(formula))
    assert interaction_signature(result.grammar) == pairs
```

- [ ] **Step 2: Run focused tests and verify failure**

Run: `pytest tests/chemistry/test_grammar.py tests/chemistry/test_analyzer.py -q`

- [ ] **Step 3: Implement frozen grammar records and stateless analyzer**

```python
@dataclass(frozen=True, slots=True)
class ChemistryResult:
    composition: Composition
    classification: ChemicalClassification
    grammar: CompositionGrammar

class ChemistryAnalyzer:
    def __init__(self, reference: ReferenceData | None = None):
        self.reference = reference or ReferenceData.default()

    def analyze(self, composition: Composition) -> ChemistryResult:
        classification = classify_composition(composition, self.reference)
        grammar = compile_composition_grammar(composition, classification, self.reference)
        return ChemistryResult(composition, classification, grammar)
```

- [ ] **Step 4: Run Chemistry and Reference Data tests and commit**

Run: `pytest tests/chemistry tests/reference_data -q`

```bash
git add src/cristma tests/chemistry
git commit -m "feat(chemistry): add interaction grammar and analyzer"
```

---

### Task 6: Package and scientific integration gate

**Files:**
- Modify: `tests/test_public_api.py`
- Create: `tests/integration/test_chemistry_core.py`
- Modify: `pyproject.toml` only if package-data discovery does not include JSON resources.

**Interfaces:**
- Verifies: `cristma.read(path) -> structure -> Composition.from_structure -> ChemistryAnalyzer.analyze`.
- Verifies wheel contains reference resources and imports without pymatgen or Qt.

- [ ] **Step 1: Add one end-to-end real-CIF test and public import test**

```python
def test_real_cif_reaches_actionable_chemistry():
    structure = cristma.read(FIXTURE).structures[0]
    result = ChemistryAnalyzer().analyze(Composition.from_structure(structure))
    assert result.classification.primary_family == "inorganic.oxide"
    assert result.grammar.candidate_interactions
```

- [ ] **Step 2: Run the new integration test**

Run: `pytest tests/integration/test_chemistry_core.py tests/test_public_api.py -q`

- [ ] **Step 3: Run the complete suite once**

Run: `pytest -q`
Expected: all tests pass.

- [ ] **Step 4: Build and inspect the wheel**

Run: `python -m build --wheel`

Inspect the wheel archive and assert both chemical reference JSON resources plus element/radius data are packaged; install it into a temporary clean virtual environment and verify imports of `cristma`, `ReferenceData`, and `ChemistryAnalyzer` with no pymatgen or Qt installed.

- [ ] **Step 5: Commit the integration gate**

```bash
git add pyproject.toml tests/integration/test_chemistry_core.py tests/test_public_api.py
git commit -m "test(chemistry): verify packaged chemistry workflow"
```
