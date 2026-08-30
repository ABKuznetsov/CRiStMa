# CRiStMa as an independent crystallographic toolbox

Date: 2026-08-30

Status: canonical architectural contract

## 1. Mission

> **CRiStMa is an independent crystallographic toolbox. It provides reusable
> scientific data types, functions, and configurable tools. Consumer context
> and orchestration remain outside its scientific API.**

CRiStMa provides one tested implementation of reusable crystallographic
mathematics for independent scientific software.

The governing boundary is:

> **Consumers own workflow and context. CRiStMa owns reusable
> crystallographic concepts and calculations.**

CRiStMa is installable and usable independently of Sci and any consumer.

## 2. Direction of dependencies

Sci is the shared platform and distribution layer for the application family.
CRiStMa is a separate scientific package distributed through that environment.

```text
Applications -> Sci
Applications -> CRiStMa
Sci          -> CRiStMa
CRiStMa      -X-> Sci or Applications
```

Sci may install and pin a tested CRiStMa release for every application. This
does not make Sci part of CRiStMa's API. CRiStMa never imports Sci, exposes Sci
types, or assumes that Sci is installed.

CRiStMa declares only dependencies required by its own scientific
implementation. Sci owns the wider application environment, common dependency
constraints, runtime validation, launchers, and deployment policy.

## 3. Package shape

The toolbox grows horizontally through independent scientific domains:

```text
cristma
|- chemistry       species, oxidation, valence knowledge
|- structure       compact canonical structure data
|- symmetry        operations, groups, orbits, Wyckoff semantics
|- io              native structure formats and optional adapters
|- geometry        distances, angles, neighbors, coordination
|- transforms      cell, setting, supercell, and entity transforms
|- hierarchy       polyhedra, units, blocks, and parameterizations
|- topology        periodic graphs, rings, dimensionality, nets
|- diffraction     reciprocal space, scattering, powder, single crystal
|- refinement      parameters, constraints, objectives, optimization helpers
`- analysis        reusable comparisons and scientific descriptors
```

The presence of one domain never makes every other domain mandatory. A powder
calculator does not need hierarchy support. A topology analyzer does not need
refinement. Consumers import only the tools they use.

Existing public namespaces such as `cristma.structure`, `cristma.symmetry`,
and `cristma.io` remain stable. Package rearrangement is not required merely to
match the conceptual diagram.

## 4. Compact scientific data types

Scientific data objects represent values and relationships. They do not become
service objects with dozens of unrelated methods.

Examples include:

```text
CrystalStructure
IndependentSite
ExpandedAtomRef
AtomicView
PeriodicNeighborGraph
PolyhedralEntity
StructuralHierarchy
ReflectionDataset
PowderPattern
CalculatedProfile
```

`CrystalStructure` remains an immutable snapshot with a cell, space group, and
independent sites. It does not acquire methods such as `refine()`,
`calculate_xrd()`, `find_topology()`, or `build_hierarchy()`.

There is no universal `ScientificObject` base class. Concrete immutable data
classes and small capability protocols are preferred over a deep inheritance
tree.

## 5. Functions and configurable tool classes

Scientific behavior lives in independent functions and tools:

```python
neighbors = NeighborFinder(method="bond_valence").find(crystal)
polyhedra = PolyhedronBuilder().build(crystal, neighbors.graph)
hierarchy = HierarchyBuilder().build(crystal)
topology = TopologyAnalyzer().analyze(crystal)
profile = PowderCalculator(radiation="CuKa").calculate(crystal, experiment)
result = RietveldRefiner(...).refine(model, observations)
```

Use a function for a simple unambiguous operation:

```python
distance = distance_between(crystal, atom_a, atom_b)
supercell_result = make_supercell(crystal, (2, 2, 1))
supercell = supercell_result.structure
```

Use a tool class when an algorithm has configuration, alternative methods, or
reusable scientific policy. Use a transform object when an operation must be
inspected, serialized, repeated, composed, inverted, or differentiated.

A tool stores inspectable configuration and remains reusable and re-entrant:

```python
finder.get_config()
finder.clone(tolerance=0.2)
```

Tool classes do not require a common superclass. Consistent conventions and
small protocols are sufficient.

## 6. Input and result contract

Inputs are treated as immutable. A call returns an explicit result instead of
mutating the input or the tool:

```text
Tool configuration + scientific input
                 |
                 v
Result
|- scientific output
|- diagnostics
`- operation mappings/provenance when scientifically required
```

Simple calculations may return a scientific value directly. Complex
calculations use a dedicated result type, for example `NeighborResult`,
`PowderCalculationResult`, or `RefinementResult`.

Diagnostics are machine-readable. Invalid API use raises an exception;
recoverable source defects, approximations, excluded observations, and
scientific warnings are represented as diagnostics.

Operation provenance records the algorithm, version, settings, assumptions,
and source-to-result mapping needed to understand the calculation.

## 7. Composition without a mandatory pipeline

Tools communicate only through typed scientific data. CRiStMa does not impose
a global pipeline such as:

```text
neighbors -> polyhedra -> hierarchy -> topology -> mechanics -> diffraction
```

The caller may choose that composition, skip steps, replace an
intermediate representation, or invoke a tool independently:

```text
                   CrystalStructure
                  /        |         \
                 v         v          v
        NeighborFinder  PowderCalculator  SymmetryAnalyzer
              |
              v
       PeriodicNeighborGraph
          /             \
         v               v
PolyhedronBuilder   TopologyAnalyzer
```

A tool may require another tool's result only when that input is a real
mathematical requirement. It depends on the result type, not on the producer
instance.

## 8. Structural transforms

Transforms are a specialized family of inspectable scientific tools:

```text
StructureTransform.apply(CrystalStructure)
                     |
                     v
TransformationResult
|- structure
|- identity_mapping
|- coordinate_mapping
|- diagnostics
`- inverse
```

Every transform preserves the origin of its result. `IdentityMapping` supports
one-to-one, one-to-many, many-to-one, removed, created, equivalent, and
ambiguous relationships. `CoordinateMapping` describes the corresponding
geometric operation.

Changing a cell must state what is preserved:

```text
preserve fractional -> fractional coordinates fixed, Cartesian geometry changes
preserve Cartesian  -> Cartesian geometry fixed, fractional coordinates change
```

Exact setting and origin transformations use exact arithmetic where possible.
Supercell and subcell operations report expanded or merged identities and
information loss explicitly.

`InverseStatus` distinguishes exact, partial, ambiguous, and unavailable
inverses. A `TransformPipeline` composes transforms and their mappings.

`DifferentiableStructureTransform` is an optional capability with evaluation
and Jacobian methods. It covers cell strain, block translation/rotation,
hinges, tilts, and internal modes used as refinement degrees of freedom.
Discrete transforms such as symmetry reduction do not pretend to be smooth
refinement parameterizations.

## 9. Hierarchy, topology, and kinematics

Hierarchy and topology are reusable scientific representations. They may
describe polyhedra, rings, blocks, chains, layers, frameworks, periodic nets,
dimensionality, hinges, and building units.

Structural entities resolve to canonical expanded atoms and independent sites.
Shared atoms remain unique. Overlap and connectivity are explicit relations,
so hierarchy is a DAG/hypergraph rather than nested copies of atoms.

Derived nodes and relations retain the mappings required to resolve them back
to the source structure. There is no universal derived-representation manager;
each tool returns the result type appropriate to its calculation.

Moving a hierarchy entity is a transform. Membership is resolved to expanded
atoms, shared-atom rules are evaluated by a kinematic model, and exactly one
new atomic coordinate is produced per atom. Incompatible controls return a
kinematic-conflict diagnostic instead of silently overwriting coordinates.

Grouping does not assert rigidity. Mechanical characterization is a separate
calculation whose result includes its assumptions and metrics. Formal scalar
rigidity scores are initially restricted to polyhedra; larger entities report
internal deformation, relative rotation, translation, and hinge changes.

## 10. Diffraction observation models

Powder and single-crystal calculations share structure-factor physics but use
different observation models:

```text
CrystalStructure + Radiation + reciprocal vectors
                       |
                       v
              StructureFactorCalculation
                       |
                       v
              ReflectionCalculation[]
                 /                 \
                v                   v
 PowderObservationModel   SingleCrystalObservationModel
                |                   |
                v                   v
 CalculatedProfile       CalculatedReflectionDataset
```

The shared layer owns reciprocal-vector generation, systematic absences,
scattering factors, occupancy, ADP, anomalous terms, and complex `F(hkl)`.
Powder multiplicity, Lorentz-polarization factors, peak overlap, and profile
convolution belong to the powder model. Measured-HKL matching, Friedel policy,
and single-crystal observation corrections belong to the single-crystal model.

A multiphase powder request may contain local `PhaseContribution` values with
a crystal and scale.

## 11. Refinement

Refinement is a reusable calculation over explicit scientific inputs:

```text
structure(s)
+ observations
+ physical model
+ parameterization
+ constraints/restraints
        |
        v
RefinementResult
|- refined scientific values
|- parameter estimates and uncertainties
|- objective statistics
|- covariance/Jacobian diagnostics
`- calculation diagnostics
```

Powder and single-crystal objectives share parameters, dependencies,
constraints, covariance tools, and optimization helpers while retaining their
own observation models.

## 12. I/O and external adapters

Native readers map CIF, RES/INS, and other structure files into compact CRiStMa
objects. They preserve source information required for scientific
round-tripping but do not discover neighboring files without an explicit
resolver.

`StructureCollection` and lazy `StructureSequence` represent finite
multi-model documents and indexed trajectory/frame sources.

Optional adapters may map to Gemmi, ASE, pymatgen, RDKit, GSAS-II, or other
ecosystems. CRiStMa objects remain the canonical inputs and outputs of its
scientific API, and every adapter remains an optional dependency.

## 13. Evolution rule

The inclusion test is:

> **Scientific + reusable + workflow-independent = CRiStMa candidate.**

A capability does not need multiple consumers before inclusion. It must be a
coherent crystallographic representation, transformation, or calculation with
an independent scientific API.

CRiStMa grows by adding new functions, compact data types, and independent
tools over the same canonical structure model. Adding `BondValenceAnalyzer`,
`TopologyAnalyzer`, `PowderCalculator`, or `SingleCrystalCalculator` does not
require redesigning existing tools or adding methods to `CrystalStructure`.

Backward-compatible public APIs are preferred. Breaking scientific semantics
requires an explicit version boundary and migration notes. Sci pins a tested
CRiStMa release so that downstream upgrades occur deliberately.

## 14. Scientific verification

Each tool is tested independently against the strongest available reference:

- analytic fixtures before external-engine comparisons;
- exact symmetry and identity mappings;
- conservation and inverse properties for transforms;
- per-reflection diagnostics for diffraction;
- numerical Jacobian checks for differentiable transforms;
- covariance and uncertainty checks for refinement;
- real structure-file fixtures with recorded provenance;
- round-trip and installed-wheel tests for public I/O.

The test suite validates CRiStMa through its public scientific API and built
distribution artifacts.

This architecture provides one implementation of crystallographic mathematics
with many independent ways to use it.
