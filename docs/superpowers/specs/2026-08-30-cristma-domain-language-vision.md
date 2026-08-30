# CRiStMa as a scientific domain language

Date: 2026-08-30

Status: architectural direction

## 1. Mission

CRiStMa is not a CIF utility, a Rietveld wrapper, or a GUI toolkit. It is the
shared scientific domain language used by independent crystallographic,
crystal-chemical, diffraction, and refinement applications.

Its purpose is to provide stable, inspectable, physics-first representations
of:

- structure and symmetry;
- measurements and experimental geometry;
- reflections, structure factors, and calculated profiles;
- parameters, dependencies, constraints, and restraints;
- derived geometry, topology, polyhedra, and structural units.

Applications consume these objects and provide workflows. They do not define
incompatible copies of the same scientific concepts.

## 2. Core boundary

The governing rule is:

> CRiStMa owns scientific entities; applications own workflows.

CRiStMa may contain:

- `CrystalStructure` and `MolecularStructure`;
- `PowderPattern` and `CalculatedProfile`;
- `Reflection` and `ReflectionSet`;
- `Parameter`, `Dependency`, `Constraint`, and `Restraint`;
- `StructureSeries` and other scientific relationships.

CRiStMa does not contain:

- project trees and recent-file state;
- Qt widgets, dialogs, or view models;
- button actions or application modes such as Guide/User/Pro;
- batch-pause and confirmation workflows;
- application-specific persistence services.

No application UI, project-state object, or service becomes part of a CRiStMa
scientific model.

## 3. Shared vocabulary

The long-lived domain vocabulary is organized into independent packages:

```text
cristma
|- structure       crystals, molecules, sites, expanded atoms
|- symmetry        operations, groups, orbits, Wyckoff semantics
|- chemistry       species, oxidation, bonds, valence knowledge
|- experiment      radiation, geometry, instruments, measurements
|- diffraction     reflections, structure factors, intensities
|- profile         observed and calculated profile components
|- refinement      parameters, dependencies, constraints, restraints
|- geometry        distances, angles, coordination, neighbor graphs
|- topology        polyhedra, structural units, blocks, hinges
`- io              native formats and optional external adapters
```

Package dependencies point toward stable domain primitives. The domain model
does not depend on readers, optimizers, application services, or rendering
frameworks.

## 4. Semantic identity

Stable semantic identifiers connect the same physical entity across readers,
viewers, diffraction, refinement, and crystal chemistry.

An independent crystallographic position contains at least:

```text
IndependentSite
|- site_id
|- label
|- species/components
|- occupancy
|- fractional_position
|- displacement
`- provenance
```

A symmetry-expanded atom is derived and contains:

```text
ExpandedAtom
|- atom_id
|- source_site_id
|- symmetry_operation_id
|- cell_translation
|- Cartesian/fractional position
`- provenance
```

Expanded atoms cannot silently become independent refinable parameters. Every
application can trace them back to the same independent site and operation.

The same principle applies to phases, reflections, measurements, profile
components, structural units, and parameters. Human-readable labels are not
used as unique identity.

## 5. Structure domain

The structure model distinguishes periodic crystallographic structures from
non-periodic molecular structures while providing a shared atomic view:

```text
Structure
|- CrystalStructure
|  |- UnitCell
|  |- SpaceGroupDefinition
|  |- IndependentSite[]
|  `- ExpandedStructure (derived)
`- MolecularStructure
   |- Atom[]
   |- Bond[]
   `- groups/residues/components
```

Derived structural representations include:

```text
CrystalStructure
-> ExpandedStructure
-> PeriodicNeighborGraph
-> CoordinationEnvironments
-> Polyhedra
-> StructuralUnits
-> RigidBlocks
-> Layers / chains / rings / hinges
```

Each derived object retains identifiers of the source atoms/sites. Structural
blocks and polyhedra therefore remain connected to refinement parameters and
visual selections.

## 6. Measurement and profile domain

Observed powder data is more than two arrays:

```text
PowderPattern
|- pattern_id
|- axis
|- observed_signal
|- sigma
|- mask
|- experiment_id
|- metadata
`- provenance
```

Experimental conditions are explicit:

```text
PowderExperiment
|- Radiation
|- MeasurementGeometry
|- InstrumentModel
|- SpecimenModel
`- PowderPattern
```

A calculated profile exposes its decomposition:

```text
CalculatedProfile
|- profile_id
|- grid
|- background
|- phase_components
|- fixed_components
|- total
|- diagnostics
`- provenance
```

Applications can inspect and visualize every component rather than receiving
only a final calculated array.

## 7. Diffraction domain

A reflection is a stable scientific entity, not a temporary row inside a
calculator:

```text
Reflection
|- reflection_id
|- phase_id
|- hkl
|- d_spacing
|- multiplicity
|- radiation_component_id
|- coordinate on the measurement axis
|- complex structure factor
|- squared amplitude
|- geometric/polarization factors
|- integrated_intensity
|- status and diagnostics
`- provenance
```

`ReflectionSet` connects a structure, experiment, and deterministic generation
settings. Finder, profile calculation, diagnostics, and refinement consume the
same reflection semantics.

## 8. Refinement domain

Refinement is expressed through scientific state rather than optimizer-specific
vectors:

```text
RefinementModel
|- structures/phases
|- experiment
|- forward model
|- Parameter[]
|- Dependency[]
|- Constraint[]
`- Restraint[]
```

A `Parameter` has a semantic owner and path, unit, value, allowed domain,
refinability, uncertainty, and provenance. Dependencies express shared or
derived parameters. Constraints are exact admissibility rules. Restraints add
scientifically justified penalty terms without pretending to be observations.

Optimizers operate through this contract. CrysPy, GSAS-II, SciPy, or a future
native optimizer may be adapters or engines; none defines CRiStMa's parameter
semantics.

Structural refinement may expose different coordinate systems over the same
underlying atoms:

```text
structural blocks -> polyhedra -> independent atoms
```

Rigid blocks and polyhedra are implemented as parameterizations and
dependencies between stable atom/site identities, not duplicated structures.

## 9. Application roles

Applications remain independent because they select only the CRiStMa domains
they need:

```text
XRD Finder
|- CrystalStructure
|- PowderPattern
`- ReflectionSet / CalculatedProfile

Rietveld Manager
|- CrystalStructure
|- PowderExperiment
|- CalculatedProfile
`- RefinementModel / RefinementState

Crystal Blocks
|- CrystalStructure
|- ExpandedStructure
|- NeighborGraph
`- StructureHierarchy
```

An application may add workflow state around these entities, but it does not
subclass or modify them with UI concerns.

## 10. Provenance and immutability

Canonical scientific objects are immutable snapshots. Transformations produce
new snapshots and record their parent state and operation. This supports:

- refinement history and undo;
- reproducible analysis;
- comparison across a temperature/pressure series;
- source-file traceability;
- safe sharing between independent applications.

Raw input, interpreted values, standard uncertainties, units, missing states,
and derived values coexist. Recovery and inference are always labeled.

## 11. Interoperability

Native format readers map external documents into CRiStMa entities. Optional
adapters map to and from Gemmi, ASE, pymatgen, RDKit, GSAS-II, and other
ecosystems.

The dependency direction is:

```text
applications -> CRiStMa
optional adapters -> CRiStMa
CRiStMa domain -X-> application UI / project state
```

An external object is never stored as the authoritative internal state. The
adapter boundary is explicit and testable.

## 12. Evolution rule

New functionality enters CRiStMa only when it can be stated as a reusable
scientific concept with stable semantics. Application-specific workflow stays
in the application.

This rule keeps CRiStMa broad in scientific capability but compact in
responsibility: one shared scientific language, many independent applications.
