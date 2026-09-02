# CRiStMa Chemistry and Reference Data Core Design

**Date:** 2026-09-01  
**Status:** approved  
**Scope:** the first reusable Chemistry and Reference Data slice for CRiStMa

## 1. Position in CRiStMa

CRiStMa is organized into seven independent scientific areas:

```text
1. I/O
2. Chemistry
3. Crystallography
4. Crystal Chemistry
5. Mathematics
6. Reference Data
7. Plugins / Workflows
```

The dependency direction for this slice is:

```text
Reference Data ──> Chemistry
       ^                |
       |                v
       +──── Crystal Chemistry <──── Crystallography

I/O ──> canonical structures
Plugins / Workflows ──> compose public tools
```

Chemistry and Reference Data do not depend on application packages, Qt,
pymatgen, CRAFT, or application workflow state.

## 2. Scientific boundary

Chemistry receives composition and immutable reference knowledge. It classifies
the chemical system and produces a composition-level grammar for later
crystallographic tests.

```text
Composition
    -> ChemicalClassification
    -> CompositionGrammar
```

Chemistry does not inspect coordinates, distances, neighbor graphs, bonds,
coordination numbers, molecules, polyhedra, structural blocks, or topology.

The boundary is:

```text
Chemistry          -> which interactions are meaningful candidates
Crystallography    -> which positions are geometrically connected
Crystal Chemistry  -> what the validated graph means chemically and structurally
```

A candidate interaction is never a confirmed bond.

## 3. Package layout

```text
cristma/
├── chemistry/
│   ├── __init__.py
│   ├── elements.py
│   ├── composition.py
│   ├── evidence.py
│   ├── classification.py
│   ├── grammar.py
│   └── analyzer.py
│
└── reference_data/
    ├── __init__.py
    ├── facade.py
    ├── elements.py
    ├── radii.py
    ├── chemical_reference.py
    └── resources/
```

Only modules required by the first implementation are created. Ionic radii,
electronegativity, oxidation-state, BVS, bond-length, and coordination datasets
do not receive placeholder modules or speculative public contracts.

## 4. Reference Data

Reference Data owns inert, versioned scientific knowledge:

```text
Chemical Reference DB
atomic properties
ionic radii
covalent radii
electronegativity scales
oxidation-state tables
BVS parameters
bond-length statistics
coordination statistics
future curated datasets
```

It does not execute classification rules, infer bonds, or analyze structures.
Algorithms belong to Chemistry or Crystal Chemistry.

Every implemented dataset exposes only:

- stable dataset and version identifiers;
- source citations and dataset provenance;
- explicit units and conventions;
- immutable records;
- exact lookup criteria;
- explicit lookup failure instead of guessed values.

If a later requirement introduces competing scientific scales, the relevant
dataset will identify its scale. The first slice does not build a general
multi-scale framework before one is needed.

### 4.1 ReferenceData facade

```python
reference = ReferenceData.default()

element = reference.elements.by_symbol("Fe")
radius = reference.covalent_radii.find("Fe")
family = reference.chemical.families["inorganic.oxide"]
```

`ReferenceData` is a small immutable bundle of these three catalogs, not a
general data platform. An analyzer may receive another compatible bundle for
tests or a later curated release. CRiStMa ships the tested default inside its
wheel and performs no network access at runtime.

### 4.2 Element records

The initial element catalog replaces Chemistry's dependency on
`pymatgen.core.Element`. At minimum an `ElementRecord` contains:

```text
symbol
atomic_number
category / metal classification
reference provenance
```

Element identity remains compatible with the existing `ElementSpecies`,
`ChargedSpecies`, `IsotopeSpecies`, and `UnknownSpecies` contracts. Display
colors are not chemical reference data and remain a Viewer policy.

### 4.3 Chemical Reference DB

The existing CRAFT v3 knowledge is migrated rather than reimplemented. The
active database is the v3.1 grammar-bearing resource. The historical v3 source
remains independently loadable for integrity and migration tests.

The supplied draft contains one invalid `CaSi2` reference:

```text
inorganic.silicide  ->  inorganic.tetrelide
```

The corrected value must refer to the existing family. Integrity validation
must reject missing family IDs, cyclic inheritance, invalid profile routes,
unknown grammar selectors or operations, malformed boundary cases, and an
unsupported schema version.

JSON remains inert data. CRiStMa never evaluates arbitrary predicates or code
stored in a reference resource.

## 5. Chemistry data model

### 5.1 Chemical entities

`chemistry/elements.py` is the base chemical-entity layer. It exposes the
existing species concepts and element normalization without importing
Reference Data or crystallographic structures.

```text
ChemicalSpecies
├── ElementSpecies
├── IsotopeSpecies
├── ChargedSpecies
└── UnknownSpecies
```

### 5.2 Composition

`Composition` is an immutable normalized chemical composition:

```text
Composition
├── amounts by normalized element symbol
├── normalized formula
└── normalization basis
```

Amounts must be positive and finite. Element symbols are normalized through
CRiStMa's own element catalog. Composition contains no coordinates or site
identity.

`Composition.from_structure(structure)` is a one-way convenience adapter over
a small structural-typing protocol. Chemistry does not import concrete
Crystallography classes, so the convenience API does not create a package
cycle. It reads only species, occupancy, and already established multiplicity:

```text
amount = occupancy * multiplicity
```

For a molecular or fully expanded atomic structure every listed position has
multiplicity one. For independent crystallographic sites, calculated
multiplicity takes precedence. An explicit identity-only structure also has
multiplicity one. A non-identity structure lacking a calculated multiplicity is
rejected rather than silently trusting a contradictory reported value.
Reported and calculated multiplicity mismatches remain visible through I/O
diagnostics.

Occupied `UnknownSpecies` cannot silently become an element. Composition
construction reports an explicit failure until the caller resolves the
identity.

### 5.3 Evidence and diagnostics

`ChemicalEvidence` records a stable code, explanatory message, elements, and
reference provenance. Chemistry reuses the toolbox-wide `Diagnostic` and
`Severity` contracts rather than defining a competing diagnostic type.

Evidence is scientific support for a hypothesis. A diagnostic describes a
problem or limitation. They are not interchangeable.

### 5.4 Classification

`ChemicalClassification` contains:

```text
primary family or unresolved
ranked alternative families
domain and composition kind
modifiers
confidence
evidence
diagnostics
reference provenance
method provenance
```

Family identities are validated reference IDs such as `inorganic.oxide`, not a
closed Python enum containing every family. Small stable concepts such as
composition kind and domain may use enums. This lets Reference Data add a
family without requiring a new Chemistry release merely to extend an enum.

The classification must be understandable at the material level. Its family
hypotheses cover at least:

```text
elemental substance
├── metallic
├── covalent
└── molecular

compound
├── intermetallic / ordered alloy candidate
├── oxide / oxysalt candidate
├── halide candidate
├── other inorganic anion family
├── organic molecular / salt candidate
├── metal-organic / organometallic candidate
└── unresolved competing families
```

The normal result contains one actionable primary family. Alternatives are
included only for an explicit curated grey-zone rule and are limited to the
few competing interpretations named by that rule. Chemistry must not return a
generic list of every theoretically possible family. For example, an elemental
substance can be classified directly as metallic, while a genuinely ambiguous
metal-carbon-hydrogen-nitrogen composition may retain organic,
metal-coordination, and organometallic alternatives until structure-derived
evidence distinguishes them.

Composition-only ambiguity is retained. Carbon, hydrogen, sulfur, phosphorus,
or a halogen does not by itself prove an organic molecule, oxoanion, molecular
anion, or confirmed mixed-anion structure.

### 5.5 Composition grammar

`CompositionGrammar` contains:

```text
decomposition mode
chemical subsystems
provisional element roles
candidate interactions
operation, layer, and priority
preferred centre/ligand direction
confidence
evidence and diagnostics
reference/method provenance
```

The initial role vocabulary includes cationic, anionic, anion subsystem,
structural former, covalent network, molecular component, metal centre,
organic component, counterion candidate, metallic, and ambiguous roles.

The initial grammar operations are the already tested composition-level
questions from CRAFT: centre-ligand shell, covalent network,
intra-subsystem bonds, interstitial coordination, mixed-anion coordination,
and metallic coordination.

Grammar records are compiled from versioned Reference Data templates. Python
does not select an interaction through material-family branches. The axes of
meaning remain independent:

```text
GrammarOperation       chemical interaction mode
InteractionLayer       role in the structure
InteractionPriority    primary or allowed chemistry request
```

Later, `ContactClassification` records primary/secondary geometric-shell
membership and must not replace any of these fields.

The practical purpose of the grammar is to tell downstream crystallographic
tools what to search for. Typical outputs are:

```text
elemental metallic -> metallic coordination candidates
intermetallic       -> heteroatomic and metallic coordination candidates,
                       without assuming cation/anion roles
oxide/oxysalt       -> structural-former-O structural candidates
                       + remaining electropositive-O interstitial candidates
halide              -> centre-halogen primary coordination candidates
organic             -> intramolecular covalent candidates first;
                       intermolecular contacts remain a separate layer
metal-organic       -> organic covalent subsystem + metal-donor candidates
```

Which concrete bonds, coordination shells, polyhedra, molecules, or networks
exist is still decided later from coordinates.

The first acceptance gate uses exact, compact outcomes rather than broad
hypothesis generation:

```text
Fe      -> elemental.metallic -> metallic coordination
Si      -> elemental.covalent -> covalent network
FeAl    -> inorganic.intermetallic -> Fe-Al/metallic coordination
CaO     -> inorganic.oxide -> Ca-O interstitial coordination candidate
NaCl    -> inorganic.halide -> Na-Cl primary coordination
FeS2    -> inorganic.chalcogenide -> Fe-S coordination + S-S intra-subsystem
```

Only curated grey-zone examples are expected to remain unresolved.

## 6. Tools and public API

Low-level functions remain independently usable:

```python
classification = classify_composition(composition, reference)
grammar = compile_composition_grammar(
    composition,
    classification,
    reference,
)
```

The convenience tool is stateless with respect to scientific inputs and
results:

```python
analyzer = ChemistryAnalyzer(reference=ReferenceData.default())
result = analyzer.analyze(composition)
```

```text
ChemistryResult
├── composition
├── classification
└── grammar
```

`ChemistryAnalyzer` stores configuration and selected references only. It has
no `current_structure`, `last_result`, or hidden application session.

## 7. Migration from CRAFT

The supplied CRAFT handoff archive (SHA-256
`d8744b29a1249237baf6cff5f8e64c65cde731b342aa037e5c405a2e60a3ce38`)
supplies migration material:

```text
chemical_reference.py      -> cristma.reference_data
chemical_classification.py -> cristma.chemistry.classification/evidence
composition_grammar.py     -> cristma.chemistry.grammar
scientific tests           -> CRiStMa tests
```

Migration is not a blind copy:

- replace `pymatgen.core.Element` with CRiStMa element/reference contracts;
- use CRiStMa `Diagnostic` and provenance types;
- separate immutable data loading from algorithms;
- remove Viewer types and application paths;
- retain grey-zone, transfer, integrity, and boundary-case tests;
- keep Stage B structural refinement in CRAFT until the future Crystal
  Chemistry design is approved.

The archive's UI, pymatgen adapter, progressive workflow, periodic-bond
resolver, hierarchy, cache, and Stage B modules are explicitly not copied into
this slice. They remain evidence for later Crystallography and Crystal
Chemistry designs.

After CRiStMa tests the extracted implementation, CRAFT imports it and deletes
its duplicate Stage A implementation. CRAFT continues to own display colors,
UI state, caching, and workflow.

## 8. First implementation slice

The first slice is deliberately complete and narrow:

1. self-contained element identity and metal classification;
2. covalent radii needed by CRAFT, with source metadata and no guessed
   silent fallback;
3. validated Chemical Reference DB v3/v3.1 resources;
4. `Composition` and its canonical-structure adapter;
5. composition classification and grammar compilation;
6. `ChemistryResult` and `ChemistryAnalyzer`;
7. migrated CRAFT Stage A scientific tests;
8. wheel installation and dependency audit proving that pymatgen and Qt are
   not required.

Ionic radii, oxidation-state solving, BVS parameters, bond-length statistics,
coordination statistics, structure-refined classification, confirmed groups,
and other Crystal Chemistry algorithms follow as independent slices over the
same contracts.

## 9. Invariants

1. Chemistry consumes composition, never coordinates.
2. Reference Data contains knowledge, never hidden algorithms.
3. Applications own workflows and display policy.
4. Every scientific result carries evidence and provenance.
5. Ambiguity is represented only for an explicit curated grey-zone rule;
   ordinary compositions produce one actionable primary family.
6. Canonical structures remain owned by Crystallography; Chemistry only derives
   a composition payload.
7. New reference datasets extend catalogs without changing Chemistry's core
   scientific objects.
