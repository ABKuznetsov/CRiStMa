# Coordination polyhedra implementation plan

1. Add contract tests for vertices, descriptors, canonical face signatures and
   degenerate geometry.
2. Extend immutable polyhedron models and builder output.
3. Build and validate space-group polyhedron orbits from oriented periodic
   vertex references.
4. Include polyhedra and their orbits in `CrystalChemistryResolution` and the
   public API.
5. Run focused, regression and package smoke verification.
