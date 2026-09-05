"""Exact symmetry-phase analysis of systematic absences."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

from cristma.crystallography import SpaceGroupSetting
from cristma.symmetry import AffineOperation

from .diagnostics import DiffractionInvariantError, INCONSISTENT_PHASE_BUCKETS
from .models import (
    ExtinctionCause,
    ExtinctionCauseKind,
    ExtinctionResult,
    MillerIndex,
    PhaseBucketEvidence,
)


FractionMatrix = tuple[tuple[Fraction, Fraction, Fraction], ...]
FractionVector = tuple[Fraction, Fraction, Fraction]


def reciprocal_action(operation: AffineOperation, hkl: MillerIndex) -> MillerIndex:
    """Apply the exact reciprocal convention ``h' = R^T h``."""

    source = hkl.as_tuple()
    values = tuple(
        sum((operation.rotation[row][column] * source[row] for row in range(3)), Fraction(0))
        for column in range(3)
    )
    if any(value.denominator != 1 for value in values):
        from .diagnostics import NON_INTEGRAL_RECIPROCAL_ACTION

        raise DiffractionInvariantError(
            NON_INTEGRAL_RECIPROCAL_ACTION,
            "symmetry rotation produced a non-integral Miller index",
            {"operation_id": operation.id, "hkl": source, "transformed": values},
        )
    return MillerIndex(*(int(value) for value in values))


def _determinant(matrix: FractionMatrix) -> Fraction:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _inverse(matrix: FractionMatrix) -> FractionMatrix:
    determinant = _determinant(matrix)
    if determinant == 0:
        raise ValueError("symmetry rotation must be invertible")
    result = []
    for row in range(3):
        values = []
        for column in range(3):
            minor_rows = [index for index in range(3) if index != column]
            minor_columns = [index for index in range(3) if index != row]
            minor = (
                matrix[minor_rows[0]][minor_columns[0]]
                * matrix[minor_rows[1]][minor_columns[1]]
                - matrix[minor_rows[0]][minor_columns[1]]
                * matrix[minor_rows[1]][minor_columns[0]]
            )
            values.append(((-1) ** (row + column)) * minor / determinant)
        result.append(tuple(values))
    return tuple(result)


def _matmul(left: FractionMatrix, right: FractionMatrix) -> FractionMatrix:
    return tuple(
        tuple(
            sum((left[row][axis] * right[axis][column] for axis in range(3)), Fraction(0))
            for column in range(3)
        )
        for row in range(3)
    )


def _matvec(matrix: FractionMatrix, vector: FractionVector) -> FractionVector:
    return tuple(
        sum((matrix[row][axis] * vector[axis] for axis in range(3)), Fraction(0))
        for row in range(3)
    )


def _rank(matrix: tuple[tuple[Fraction, ...], ...]) -> int:
    rows = [list(row) for row in matrix]
    rank = 0
    row_count = len(rows)
    column_count = len(rows[0]) if rows else 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(rank, row_count) if rows[row][column] != 0),
            None,
        )
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        divisor = rows[rank][column]
        rows[rank] = [value / divisor for value in rows[rank]]
        for row in range(row_count):
            if row == rank:
                continue
            factor = rows[row][column]
            rows[row] = [a - factor * b for a, b in zip(rows[row], rows[rank], strict=True)]
        rank += 1
    return rank


def _relative_operation(reference: AffineOperation, operation: AffineOperation) -> tuple[FractionMatrix, FractionVector]:
    inverse = _inverse(reference.rotation)
    rotation = _matmul(inverse, operation.rotation)
    translation = _matvec(
        inverse,
        tuple(right - left for left, right in zip(reference.translation, operation.translation, strict=True)),
    )
    return rotation, tuple(value % 1 for value in translation)


def _has_intrinsic_translation(rotation: FractionMatrix, translation: FractionVector) -> bool:
    origin_shift_matrix = tuple(
        tuple(Fraction(int(row == column)) - rotation[row][column] for column in range(3))
        for row in range(3)
    )
    augmented = tuple(
        row + (translation[index],)
        for index, row in enumerate(origin_shift_matrix)
    )
    return _rank(augmented) > _rank(origin_shift_matrix)


@lru_cache(maxsize=None)
def _classify_pair(reference: AffineOperation, operation: AffineOperation) -> ExtinctionCauseKind:
    rotation, translation = _relative_operation(reference, operation)
    identity = tuple(
        tuple(Fraction(int(row == column)) for column in range(3))
        for row in range(3)
    )
    non_lattice = any(value % 1 for value in translation)
    if rotation == identity and non_lattice:
        return ExtinctionCauseKind.CENTERING
    nullity = 3 - _rank(
        tuple(
            tuple(rotation[row][column] - identity[row][column] for column in range(3))
            for row in range(3)
        )
    )
    intrinsic = _has_intrinsic_translation(rotation, translation)
    determinant = _determinant(rotation)
    if determinant == 1 and nullity == 1 and intrinsic:
        return ExtinctionCauseKind.SCREW_AXIS
    if determinant == -1 and nullity == 2 and intrinsic:
        return ExtinctionCauseKind.GLIDE_PLANE
    return ExtinctionCauseKind.COMBINED


@dataclass(frozen=True, slots=True)
class ExtinctionAnalyzer:
    """Determine systematic absence from exact space-group operations."""

    def analyze(self, hkl: MillerIndex, space_group: SpaceGroupSetting) -> ExtinctionResult:
        if not isinstance(hkl, MillerIndex):
            raise TypeError("hkl must be MillerIndex")
        if hkl.is_zero:
            raise ValueError("zero Miller index is not a reflection")
        if not isinstance(space_group, SpaceGroupSetting):
            raise TypeError("space_group must be SpaceGroupSetting")
        operation_ids = tuple(operation.id for operation in space_group.symmetry_operations)
        if any(not isinstance(value, str) or not value.strip() for value in operation_ids):
            raise ValueError("every symmetry operation must have a non-empty ID")
        if len(set(operation_ids)) != len(operation_ids):
            raise ValueError("symmetry operation IDs must be unique within a setting")

        grouped: dict[MillerIndex, list[tuple[AffineOperation, Fraction]]] = defaultdict(list)
        for operation in sorted(space_group.symmetry_operations, key=lambda item: item.id or ""):
            transformed = reciprocal_action(operation, hkl)
            phase = sum(
                (Fraction(index) * translation for index, translation in zip(hkl.as_tuple(), operation.translation, strict=True)),
                Fraction(0),
            ) % 1
            grouped[transformed].append((operation, phase))

        evidence = []
        bucket_operations: dict[MillerIndex, tuple[AffineOperation, ...]] = {}
        for transformed in sorted(grouped):
            items = tuple(sorted(grouped[transformed], key=lambda item: item[0].id or ""))
            reference_phase = items[0][1]
            relative = tuple((phase - reference_phase) % 1 for _operation, phase in items)
            operations = tuple(operation for operation, _phase in items)
            bucket_operations[transformed] = operations
            evidence.append(
                PhaseBucketEvidence(
                    transformed_hkl=transformed,
                    operation_ids=tuple(operation.id or "" for operation in operations),
                    translation_parts=tuple(
                        tuple(value % 1 for value in operation.translation) for operation in operations
                    ),
                    exact_phases=tuple(phase for _operation, phase in items),
                    relative_phases=relative,
                    cancels=any(phase != 0 for phase in relative),
                )
            )

        verdicts = {item.cancels for item in evidence}
        if len(verdicts) != 1:
            raise DiffractionInvariantError(
                INCONSISTENT_PHASE_BUCKETS,
                "symmetry phase buckets give inconsistent extinction verdicts",
                {
                    "hkl": hkl.as_tuple(),
                    "space_group_setting_id": space_group.setting_id,
                    "buckets": tuple(
                        {
                            "transformed_hkl": item.transformed_hkl.as_tuple(),
                            "operation_ids": item.operation_ids,
                            "exact_phases": item.exact_phases,
                            "verdict": item.cancels,
                        }
                        for item in evidence
                    ),
                },
            )

        absent = verdicts.pop()
        exact_evidence = tuple(evidence)
        if not absent:
            return ExtinctionResult(False, (), exact_evidence)

        kinds: set[ExtinctionCauseKind] = set()
        involved: set[str] = set()
        cancelled_evidence = tuple(item for item in exact_evidence if item.cancels)
        for item in cancelled_evidence:
            operations = bucket_operations[item.transformed_hkl]
            reference = operations[0]
            involved.update(item.operation_ids)
            for operation, phase in zip(operations[1:], item.relative_phases[1:], strict=True):
                if phase != 0:
                    kinds.add(_classify_pair(reference, operation))
        kind = next(iter(kinds)) if len(kinds) == 1 else ExtinctionCauseKind.COMBINED
        phase_text = ",".join(
            str(phase) for item in cancelled_evidence for phase in item.relative_phases if phase != 0
        )
        condition = f"h={hkl.as_tuple()}: non-trivial exact phase character ({phase_text})"
        cause = ExtinctionCause(kind, tuple(sorted(involved)), cancelled_evidence, condition)
        return ExtinctionResult(True, (cause,), exact_evidence)


__all__ = ["ExtinctionAnalyzer"]
