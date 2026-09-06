from __future__ import annotations

from fractions import Fraction
from itertools import combinations_with_replacement, product
import math

import numpy as np

from cristma.core import MeasuredValue, UnitCell
from cristma.crystallography import (
    AsymmetricUnitMapper,
    SymmetryContext,
    SymmetryPairFinder,
)
from cristma.structure import CrystalStructure, IndependentSite, SiteComponent
from cristma.symmetry import AffineOperation


def _value(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def _cell(
    a: float,
    b: float,
    c: float,
    alpha: float = 90.0,
    beta: float = 90.0,
    gamma: float = 90.0,
) -> UnitCell:
    return UnitCell(
        _value(a),
        _value(b),
        _value(c),
        _value(alpha),
        _value(beta),
        _value(gamma),
    )


def _site(site_id: str, fractional: tuple[float, float, float]) -> IndependentSite:
    return IndependentSite(
        id=site_id,
        label=site_id,
        components=(SiteComponent("C", _value(1.0)),),
        fractional=tuple(_value(value) for value in fractional),
    )


def _fixture(cell: UnitCell, *sites: IndependentSite):
    identity = AffineOperation(
        (
            (Fraction(1), Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(1), Fraction(0)),
            (Fraction(0), Fraction(0), Fraction(1)),
        ),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    structure = CrystalStructure("pair-search", cell, tuple(sites))
    context = SymmetryContext.from_operations((identity,), cell)
    mapping = AsymmetricUnitMapper().build(structure, context)
    return structure, context, mapping


def _candidate_descriptor(candidate) -> tuple[str, str, tuple[int, int, int], float]:
    return (
        candidate.first_site_id,
        candidate.second_site_id,
        candidate.relation.lattice_translation,
        round(candidate.distance, 12),
    )


def _brute_force_descriptors(
    structure: CrystalStructure,
    cutoff: float,
    translation_limit: int,
) -> set[tuple[str, str, tuple[int, int, int], float]]:
    rows = set()
    by_id = sorted(structure.sites, key=lambda site: site.id)
    for first, second in combinations_with_replacement(by_id, 2):
        first_fractional = np.asarray(tuple(value.value for value in first.fractional), dtype=float)
        second_fractional = np.asarray(tuple(value.value for value in second.fractional), dtype=float)
        for translation in product(
            range(-translation_limit, translation_limit + 1),
            repeat=3,
        ):
            vector = (second_fractional + np.asarray(translation) - first_fractional) @ structure.cell.matrix
            distance = float(np.linalg.norm(vector))
            if 1e-12 < distance <= cutoff + 1e-12:
                rows.add((first.id, second.id, translation, round(distance, 12)))
    return rows


def test_skew_cell_buffer_matches_large_bruteforce_reference() -> None:
    cell = _cell(2.0, 2.2, 2.4, 70.0, 80.0, 75.0)
    structure, context, mapping = _fixture(
        cell,
        _site("A", (0.08, 0.17, 0.29)),
        _site("B", (0.81, 0.64, 0.53)),
    )

    found = SymmetryPairFinder(cutoff=2.1).find_candidates(structure, context, mapping)

    assert found.complete
    assert {_candidate_descriptor(item) for item in found.candidates} == _brute_force_descriptors(
        structure,
        2.1,
        5,
    )


def test_cutoff_boundary_uses_documented_metric_tolerance() -> None:
    cell = _cell(4.0, 4.0, 4.0)
    structure, context, mapping = _fixture(
        cell,
        _site("A", (0.0, 0.0, 0.0)),
        _site("B", (0.5, 0.0, 0.0)),
    )

    found = SymmetryPairFinder(cutoff=2.0, distance_tolerance=1e-12).find_candidates(
        structure,
        context,
        mapping,
    )

    boundary = tuple(
        item for item in found.candidates if math.isclose(item.distance, 2.0, abs_tol=1e-12)
    )
    assert len(boundary) == 2
    assert {item.relation.lattice_translation for item in boundary} == {(-1, 0, 0), (0, 0, 0)}


def test_complete_search_can_require_translation_outside_one_cell() -> None:
    cell = _cell(1.0, 1.0, 1.0)
    structure, context, mapping = _fixture(cell, _site("A", (0.0, 0.0, 0.0)))

    found = SymmetryPairFinder(cutoff=2.01).find_candidates(structure, context, mapping)

    translations = {item.relation.lattice_translation for item in found.candidates}
    assert (2, 0, 0) in translations
    assert (-2, 0, 0) in translations


def test_max_candidates_returns_controlled_incomplete_result() -> None:
    cell = _cell(1.0, 1.0, 1.0)
    structure, context, mapping = _fixture(cell, _site("A", (0.0, 0.0, 0.0)))

    found = SymmetryPairFinder(cutoff=2.01, max_candidates=3).find_candidates(
        structure,
        context,
        mapping,
    )

    assert not found.complete
    assert any(
        diagnostic.code == "symmetry.pairs.search_limit_reached"
        for diagnostic in found.diagnostics
    )
    assert dict(found.provenance)["candidate_pairs_tested"] == 3
