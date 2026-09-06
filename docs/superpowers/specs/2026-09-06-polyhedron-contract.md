# Coordination polyhedra contract

CrIStMa computes immutable scientific geometry. Presentation, mesh
triangulation, colour, visibility, selection and motif comparison belong to
consumers such as CRAFT.

## Result

`CrystalChemistryResolution` exposes individual coordination polyhedra and
their space-group orbits. Each polyhedron contains stable identities, oriented
periodic vertex references, coordination number, effective ligand composition,
convex-hull faces, metric descriptors, scientific status, diagnostics and
provenance.

`polyhedron_orbit_id` is derived from the action of the supplied space-group
operations on the centre and every oriented periodic vertex reference. It is
not inferred from approximate shape similarity.

## Descriptors

For bond lengths `d_i`, the Baur distortion index is

`D = mean(abs(d_i - mean(d)) / mean(d))`.

Hull edges are the unique unordered vertex pairs occurring on face boundaries.
`edge_angle_dispersion_deg` is the population standard deviation (`ddof=0`) of
the centre-ligand-centre angles for those edges. It is comparable only when
coordination number and canonical face signature match.

The face signature canonically represents the complete vertex-edge-face
incidence graph, independent of vertex numbering. It is not merely a sorted
list of face sizes.

Degenerate or incomplete geometry is never guessed. An unavailable descriptor
is `None` and its reason is recorded in diagnostics.

## Occupation

A vertex occupancy is the effective total occupation of that geometric ligand
position. The complete component composition remains reachable by resolving
the vertex's `PeriodicAtomRef` against the atomic view. Ligand composition is
the occupancy-weighted element summary over vertices.
