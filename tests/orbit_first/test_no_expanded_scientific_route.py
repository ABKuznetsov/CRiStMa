from __future__ import annotations

import inspect
import re

from cristma.crystal_chemistry import (
    PeriodicConnectivityAnalyzer,
    PolyhedronOrbitBuilder,
    RingFinder,
    StructuralGraphBuilder,
)


def test_orbit_first_scientific_entry_points_do_not_accept_materialized_contacts() -> None:
    for method in (
        PolyhedronOrbitBuilder.build,
        StructuralGraphBuilder.build,
        PeriodicConnectivityAnalyzer.analyze,
        RingFinder.find,
    ):
        parameters = set(inspect.signature(method).parameters)
        assert "contacts" not in parameters
        assert "structure" not in parameters
        assert "atomic_view" not in parameters
        assert not any(name.startswith("legacy_") for name in parameters)


def test_scientific_modules_do_not_import_materialized_contact_type() -> None:
    modules = (
        "orbit_contacts", "incidence_orbits", "shell_orbits", "polyhedron_orbits",
        "structural_units", "structural_graph", "periodic_connectivity", "ring_finder",
    )
    for name in modules:
        module = __import__("cristma.crystal_chemistry." + name, fromlist=(name,))
        assert re.search(r"\bResolvedContact\b", inspect.getsource(module)) is None
