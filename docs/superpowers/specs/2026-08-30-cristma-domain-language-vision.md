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

CRiStMa is distributed as a standalone Python package through PyPI and is
installable with the standard Python packaging toolchain:

```bash
pip install cristma
```

## 2. Distribution and direction of dependencies

CRiStMa has a public, consumer-neutral Python API. Any script, notebook,
scientific package, desktop application, or service may depend on it directly.

```text
Python consumer -> CRiStMa
CRiStMa         -X-> consumer
```

Published releases use semantic versions and declare only dependencies required
by CRiStMa's own scientific implementation. Large integrations and specialized
format bridges are provided through optional extras or separate adapter
packages, so the base installation remains usable on its own.

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

Domains are not assigned to particular applications. Finder is not restricted
to diffraction, CRAFT is not restricted to geometry, and Rietveld software is
not restricted to refinement. Any consumer may compose any public CRiStMa
function or tool whose scientific contract fits its task. Adding a new public
tool makes it available to every consumer without an application-specific
registration step or a library change for that consumer.

> **CRiStMa does not distribute capabilities among applications. It provides a
> common scientific toolbox; each consumer selects and orchestrates the tools
> it needs.**

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
ExpandedAtom
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
view = expand_structure(crystal)
neighbors = NeighborFinder(method="bond_valence").find(view)
polyhedra = PolyhedronBuilder().build(view, neighbors)
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

A tool never stores the current structure, last result, or an implicit
scientific session:

```text
caller-owned input + tool configuration -> explicit result
```

The caller owns calculation order, reuse, invalidation, and caching of
`AtomicView`, neighbor graphs, polyhedra, profiles, and other results. Reusing a
previous `PeriodicNeighborGraph` is an explicit argument passed to the next
tool, not hidden state shared between tool instances.

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

The forbidden session-style API is:

```python
cristma.load("sample.cif")
cristma.current_structure
cristma.find_neighbors()
```

The supported style keeps ownership visible:

```python
crystal = cristma.read("sample.cif").structures[0]
view = expand_structure(crystal)
graph = NeighborFinder(cutoff=3.0).find(view)
```

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

### Canonical scientific input invariant

`CrystalStructure` and `MolecularStructure` are canonical internal scientific
models, not canonical file formats. Every format-specific representation ends
at the I/O mapper:

```text
CIF / RES / INS / POSCAR / QE / other source
        -> format document
        -> semantic mapper
        -> CrystalStructure | MolecularStructure
```

After successful mapping, source format does not select or alter scientific
mathematics. Two canonical objects representing the same structure must produce
equivalent results from symmetry, geometry, bond-valence, diffraction, and
other structure-based tools regardless of whether they originated in CIF,
SHELX, POSCAR, or another source.

Derived representations retain provenance back to the canonical structure but
never become an alternative source of truth:

```text
CrystalStructure
        -> AtomicView
        -> PeriodicNeighborGraph
        -> PolyhedronSet / StructuralHierarchy / TopologyResult
```

Refinement has no separate Rietveld-specific structure. A parameterization
produces a new `CrystalStructure`, and the same forward calculators evaluate
that structure.

> **`CrystalStructure` / `MolecularStructure` are the canonical scientific
> inputs for all structure-based calculations. File-specific representations
> terminate at the I/O boundary.**

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
requires an explicit version boundary and migration notes. Release notes state
scientific as well as API changes so consumers can pin or upgrade deliberately.

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
