from __future__ import annotations

from dataclasses import replace
from fractions import Fraction

from cristma.chemistry import InteractionLayer
from cristma.core import MeasuredValue, UnitCell
from cristma.crystal_chemistry import (
    PeriodicConnectivityAnalyzer,
    RingFinder,
    ShellRole,
    StructuralBlockFinder,
    StructuralConnectionKind,
    StructuralConnectionOrbit,
    StructuralRepresentation,
    StructuralSelectionPolicy,
    StructuralUnitKind,
    StructuralUnitOrbit,
)
from cristma.crystallography import PeriodicSymmetryRelation, SymmetryContext
from cristma.structure import PeriodicAtomRef
from cristma.symmetry import AffineOperation


IDENTITY = AffineOperation(
    ((Fraction(1), Fraction(0), Fraction(0)),
     (Fraction(0), Fraction(1), Fraction(0)),
     (Fraction(0), Fraction(0), Fraction(1))),
    (Fraction(0), Fraction(0), Fraction(0)),
)


def _value(value: float) -> MeasuredValue:
    return MeasuredValue(value, None, str(value))


def _triangle(kind: StructuralConnectionKind = StructuralConnectionKind.SHARED_VERTEX):
    cell = UnitCell(*(_value(10.0) for _ in range(3)), *(_value(90.0) for _ in range(3)))
    context = SymmetryContext.from_operations((IDENTITY,), cell)
    relation = PeriodicSymmetryRelation(context.identity_operation_key, (0, 0, 0))
    units = tuple(
        StructuralUnitOrbit(
            name, StructuralUnitKind.ATOM, name, None, None,
            (PeriodicAtomRef("image:" + name, (0, 0, 0)),), (), 1,
            (("fixture", "triangle"),),
        )
        for name in ("A", "B", "C")
    )
    roles = (ShellRole.PRIMARY,) if kind is not StructuralConnectionKind.DIRECT_CONTACT else ()
    edges = tuple(
        StructuralConnectionOrbit(
            "edge:" + first + second, first, second, relation, kind,
            ((first + second),) if kind is not StructuralConnectionKind.DIRECT_CONTACT else (),
            (InteractionLayer.STRUCTURAL,), roles, ("contact:" + first + second,),
            ("interpretation:" + first + second,), 1,
        )
        for first, second in (("A", "B"), ("B", "C"), ("A", "C"))
    )
    policy = StructuralSelectionPolicy(
        frozenset({InteractionLayer.STRUCTURAL}), frozenset({ShellRole.PRIMARY})
    )
    representation = StructuralRepresentation(
        "triangle", units, edges, policy, _symmetry_context=context
    )
    connectivity = PeriodicConnectivityAnalyzer().analyze(representation)
    blocks = StructuralBlockFinder().find(representation, connectivity)
    return representation, blocks


def test_ring_finder_returns_one_quotient_ring_orbit() -> None:
    representation, blocks = _triangle()

    result = RingFinder().find(representation, blocks)

    assert len(result.ring_orbits) == 1
    ring = result.ring_orbits[0]
    assert ring.size == 3
    assert ring.multiplicity_in_reference_cell == 1
    assert not hasattr(result, "rings")
    assert all(ref.periodic_relation.operation_key for ref in ring.unit_orbit_refs)


def test_direct_contacts_do_not_create_structural_rings() -> None:
    representation, blocks = _triangle(StructuralConnectionKind.DIRECT_CONTACT)

    assert RingFinder().find(representation, blocks).ring_orbits == ()


def test_ring_orbit_identity_is_independent_of_operation_order() -> None:
    representation, blocks = _triangle()
    inversion = AffineOperation(
        ((Fraction(-1), Fraction(0), Fraction(0)),
         (Fraction(0), Fraction(-1), Fraction(0)),
         (Fraction(0), Fraction(0), Fraction(-1))),
        (Fraction(0), Fraction(0), Fraction(0)),
    )
    cell = UnitCell(*(_value(10.0) for _ in range(3)), *(_value(90.0) for _ in range(3)))
    first = SymmetryContext.from_operations((IDENTITY, inversion), cell)
    second = SymmetryContext.from_operations((inversion, IDENTITY), cell)

    first_result = RingFinder().find(replace(representation, _symmetry_context=first), blocks)
    second_result = RingFinder().find(replace(representation, _symmetry_context=second), blocks)

    assert first_result.ring_orbits == second_result.ring_orbits
    assert first_result.ring_orbits[0].multiplicity_in_reference_cell == 2
