"""Validated, deterministic direct-space symmetry contexts."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import StrEnum
from fractions import Fraction
import hashlib
import json
import math
from types import MappingProxyType
from typing import Iterable, Mapping

import numpy as np

from cristma.core import UnitCell
from cristma.diagnostics import Diagnostic, Severity
from cristma.symmetry.affine import AffineOperation, Matrix3, Vector3
from cristma.symmetry.orbit import SpaceGroupDefinition

from .space_group import SpaceGroupSetting


DEFAULT_METRIC_TOLERANCE = 1e-8

_IDENTITY_ROTATION: Matrix3 = (
    (Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(1), Fraction(0)),
    (Fraction(0), Fraction(0), Fraction(1)),
)
_ZERO_TRANSLATION: Vector3 = (Fraction(0), Fraction(0), Fraction(0))


class SymmetrySourceKind(StrEnum):
    """How a validated direct-space symmetry action was supplied."""

    CATALOG_SETTING = "catalog_setting"
    VALID_EXPLICIT_OPERATIONS = "valid_explicit_operations"
    EXPLICIT_IDENTITY_FALLBACK = "explicit_identity_fallback"


class DirectBasisConvention(StrEnum):
    """Coordinate convention used by the direct-space affine action."""

    FRACTIONAL_DIRECT = "fractional_direct"


class SymmetryContextInvariantError(ValueError):
    """A raw symmetry source cannot form a valid working context."""

    def __init__(
        self,
        code: str,
        message: str,
        *,
        evidence: tuple[tuple[str, object], ...] = (),
    ) -> None:
        super().__init__(message)
        self.code = code
        self.evidence = evidence


def _fraction_descriptor(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _operation_descriptor(operation: AffineOperation) -> str:
    normalized = operation.normalized()
    values = (
        *(value for row in normalized.rotation for value in row),
        *normalized.translation,
    )
    return ",".join(_fraction_descriptor(value) for value in values)


def canonical_operation_key(operation: AffineOperation) -> str:
    """Return an exact key independent of source text and input ordering."""

    if not isinstance(operation, AffineOperation):
        raise TypeError("operation must be AffineOperation")
    descriptor = _operation_descriptor(operation).encode("ascii")
    return "operation:" + hashlib.sha256(descriptor).hexdigest()


def _canonical_operation(operation: AffineOperation) -> AffineOperation:
    if not isinstance(operation, AffineOperation):
        raise TypeError("operations must contain only AffineOperation values")
    rotation: list[tuple[Fraction, Fraction, Fraction]] = []
    for row in operation.rotation:
        if len(row) != 3:
            raise SymmetryContextInvariantError(
                "symmetry.context.group_invalid",
                "symmetry rotations must be 3x3 matrices",
            )
        converted = tuple(Fraction(value) for value in row)
        if any(value.denominator != 1 for value in converted):
            raise SymmetryContextInvariantError(
                "symmetry.context.group_invalid",
                "symmetry rotations must contain exact integers",
            )
        rotation.append(converted)
    if len(rotation) != 3 or len(operation.translation) != 3:
        raise SymmetryContextInvariantError(
            "symmetry.context.group_invalid",
            "symmetry operations must act in three dimensions",
        )
    translation = tuple(Fraction(value) % 1 for value in operation.translation)
    return AffineOperation(tuple(rotation), translation)


def _determinant(matrix: Matrix3) -> Fraction:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def _matrix_product(left: Matrix3, right: Matrix3) -> Matrix3:
    return tuple(
        tuple(
            sum((left[row][inner] * right[inner][column] for inner in range(3)), Fraction(0))
            for column in range(3)
        )
        for row in range(3)
    )


def _matrix_vector(matrix: Matrix3, vector: Vector3) -> Vector3:
    return tuple(
        sum((matrix[row][column] * vector[column] for column in range(3)), Fraction(0))
        for row in range(3)
    )


def _compose(left: AffineOperation, right: AffineOperation) -> AffineOperation:
    rotated_translation = _matrix_vector(left.rotation, right.translation)
    return AffineOperation(
        _matrix_product(left.rotation, right.rotation),
        tuple(
            (rotated_translation[index] + left.translation[index]) % 1
            for index in range(3)
        ),
    )


def _cell_fingerprint(cell: UnitCell) -> str:
    measured = (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    values = tuple(float(item.value) for item in measured)
    return hashlib.sha256("|".join(value.hex() for value in values).encode("ascii")).hexdigest()


def _stable_value(value: object) -> object:
    if value is None or isinstance(value, (bool, int, str)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("symmetry provenance floats must be finite")
        return {"float": value.hex()}
    if isinstance(value, Fraction):
        return {"fraction": [value.numerator, value.denominator]}
    if isinstance(value, StrEnum):
        return value.value
    if isinstance(value, Mapping):
        return {
            str(key): _stable_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_stable_value(item) for item in value]
    raise TypeError(f"unsupported symmetry provenance value: {type(value).__name__}")


def _digest(value: object) -> str:
    serialized = json.dumps(
        _stable_value(value),
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")
    return hashlib.sha256(serialized).hexdigest()


def _validate_metric(
    operations: tuple[AffineOperation, ...],
    cell: UnitCell,
    tolerance: float,
) -> None:
    metric = np.asarray(cell.metric, dtype=float)
    scale = max(1.0, float(np.max(np.abs(metric))))
    for operation in operations:
        rotation = np.asarray(operation.rotation, dtype=float)
        transformed = rotation.T @ metric @ rotation
        if not np.allclose(
            transformed,
            metric,
            rtol=tolerance,
            atol=tolerance * scale,
        ):
            raise SymmetryContextInvariantError(
                "symmetry.context.metric_incompatible",
                "symmetry rotation is incompatible with the supplied cell metric",
                evidence=(("operation_key", canonical_operation_key(operation)),),
            )


def _validate_and_canonicalize_operations(
    operations: Iterable[AffineOperation],
    cell: UnitCell,
    metric_tolerance: float,
) -> tuple[AffineOperation, ...]:
    if not isinstance(cell, UnitCell):
        raise TypeError("cell must be UnitCell")
    if isinstance(metric_tolerance, bool) or not isinstance(metric_tolerance, (int, float)):
        raise TypeError("metric_tolerance must be a real number")
    metric_tolerance = float(metric_tolerance)
    if not math.isfinite(metric_tolerance) or metric_tolerance <= 0:
        raise ValueError("metric_tolerance must be positive and finite")

    canonical = tuple(_canonical_operation(operation) for operation in operations)
    if not canonical:
        raise SymmetryContextInvariantError(
            "symmetry.context.group_invalid",
            "symmetry context requires at least one operation",
        )
    descriptors = tuple(_operation_descriptor(operation) for operation in canonical)
    if len(set(descriptors)) != len(descriptors):
        raise SymmetryContextInvariantError(
            "symmetry.context.duplicate_operation",
            "symmetry context contains duplicate normalized operations",
        )
    ordered = tuple(
        operation
        for _, operation in sorted(zip(descriptors, canonical, strict=True), key=lambda item: item[0])
    )
    for operation in ordered:
        if _determinant(operation.rotation) not in {Fraction(-1), Fraction(1)}:
            raise SymmetryContextInvariantError(
                "symmetry.context.group_invalid",
                "symmetry rotations must have determinant +1 or -1",
                evidence=(("operation_key", canonical_operation_key(operation)),),
            )

    descriptor_set = {_operation_descriptor(operation) for operation in ordered}
    identity = AffineOperation(_IDENTITY_ROTATION, _ZERO_TRANSLATION)
    identity_descriptor = _operation_descriptor(identity)
    if identity_descriptor not in descriptor_set:
        raise SymmetryContextInvariantError(
            "symmetry.context.group_invalid",
            "symmetry operation set has no identity",
        )
    for left in ordered:
        for right in ordered:
            product_descriptor = _operation_descriptor(_compose(left, right))
            if product_descriptor not in descriptor_set:
                raise SymmetryContextInvariantError(
                    "symmetry.context.group_invalid",
                    "symmetry operation set is not closed",
                    evidence=(
                        ("left_operation_key", canonical_operation_key(left)),
                        ("right_operation_key", canonical_operation_key(right)),
                    ),
                )
    for operation in ordered:
        if not any(
            _operation_descriptor(_compose(operation, candidate)) == identity_descriptor
            and _operation_descriptor(_compose(candidate, operation)) == identity_descriptor
            for candidate in ordered
        ):
            raise SymmetryContextInvariantError(
                "symmetry.context.group_invalid",
                "symmetry operation has no two-sided inverse",
                evidence=(("operation_key", canonical_operation_key(operation)),),
            )

    _validate_metric(ordered, cell, metric_tolerance)
    return ordered


@dataclass(frozen=True, slots=True)
class SymmetryContext:
    """A validated exact symmetry action in one numerical direct basis."""

    operations: tuple[AffineOperation, ...]
    operation_keys: tuple[str, ...]
    basis_convention: DirectBasisConvention
    cell_fingerprint: str
    symmetry_action_fingerprint: str
    fingerprint: str
    setting_id: str | None
    source_kind: SymmetrySourceKind
    metric_tolerance: float
    diagnostics: tuple[Diagnostic, ...]
    provenance: tuple[tuple[str, object], ...]
    _operation_lookup: Mapping[str, AffineOperation] = field(
        repr=False,
        compare=False,
        hash=False,
    )

    @classmethod
    def _build(
        cls,
        operations: Iterable[AffineOperation],
        cell: UnitCell,
        *,
        setting_id: str | None,
        source_kind: SymmetrySourceKind,
        diagnostics: tuple[Diagnostic, ...],
        provenance: tuple[tuple[str, object], ...],
        metric_tolerance: float,
    ) -> "SymmetryContext":
        canonical = _validate_and_canonicalize_operations(
            operations,
            cell,
            metric_tolerance,
        )
        keys = tuple(canonical_operation_key(operation) for operation in canonical)
        basis = DirectBasisConvention.FRACTIONAL_DIRECT
        action_fingerprint = _digest(
            {
                "basis_convention": basis.value,
                "operations": tuple(_operation_descriptor(operation) for operation in canonical),
            }
        )
        cell_digest = _cell_fingerprint(cell)
        fingerprint = _digest(
            {
                "symmetry_action_fingerprint": action_fingerprint,
                "cell_fingerprint": cell_digest,
                "setting_id": setting_id,
                "source_kind": source_kind.value,
                "metric_tolerance": metric_tolerance,
                "diagnostics": tuple(
                    (item.severity.value, item.code, item.message, item.recovery)
                    for item in diagnostics
                ),
                "provenance": provenance,
            }
        )
        return cls(
            operations=canonical,
            operation_keys=keys,
            basis_convention=basis,
            cell_fingerprint=cell_digest,
            symmetry_action_fingerprint=action_fingerprint,
            fingerprint=fingerprint,
            setting_id=setting_id,
            source_kind=source_kind,
            metric_tolerance=metric_tolerance,
            diagnostics=diagnostics,
            provenance=provenance,
            _operation_lookup=MappingProxyType(dict(zip(keys, canonical, strict=True))),
        )

    @classmethod
    def from_operations(
        cls,
        operations: Iterable[AffineOperation],
        cell: UnitCell,
        *,
        provenance: tuple[tuple[str, object], ...] = (),
        metric_tolerance: float = DEFAULT_METRIC_TOLERANCE,
    ) -> "SymmetryContext":
        diagnostics = (
            Diagnostic(
                Severity.INFO,
                "symmetry.setting_unresolved",
                "valid explicit symmetry operations are not tied to a catalog setting",
            ),
        )
        return cls._build(
            operations,
            cell,
            setting_id=None,
            source_kind=SymmetrySourceKind.VALID_EXPLICIT_OPERATIONS,
            diagnostics=diagnostics,
            provenance=tuple(provenance),
            metric_tolerance=metric_tolerance,
        )

    @classmethod
    def from_definition(
        cls,
        definition: SpaceGroupDefinition,
        cell: UnitCell,
        *,
        metric_tolerance: float = DEFAULT_METRIC_TOLERANCE,
    ) -> "SymmetryContext":
        if not isinstance(definition, SpaceGroupDefinition):
            raise TypeError("definition must be SpaceGroupDefinition")
        fallback = definition.provenance in {"identity_fallback", "unreported_identity"}
        source_kind = (
            SymmetrySourceKind.EXPLICIT_IDENTITY_FALLBACK
            if fallback
            else SymmetrySourceKind.VALID_EXPLICIT_OPERATIONS
        )
        diagnostics = (
            Diagnostic(
                Severity.WARNING if fallback else Severity.INFO,
                "symmetry.identity_fallback" if fallback else "symmetry.setting_unresolved",
                (
                    "an explicitly requested identity symmetry fallback is in use"
                    if fallback
                    else "valid explicit symmetry operations are not tied to a catalog setting"
                ),
            ),
        )
        provenance: tuple[tuple[str, object], ...] = (
            ("source", "space_group_definition"),
            ("reported_provenance", definition.provenance),
            ("number", definition.number),
            ("hall_symbol", definition.hall_symbol),
            ("setting", definition.setting),
            ("origin_choice", definition.origin_choice),
        )
        return cls._build(
            definition.operations,
            cell,
            setting_id=None,
            source_kind=source_kind,
            diagnostics=diagnostics,
            provenance=provenance,
            metric_tolerance=metric_tolerance,
        )

    @classmethod
    def from_setting(
        cls,
        setting: SpaceGroupSetting,
        cell: UnitCell,
        *,
        metric_tolerance: float = DEFAULT_METRIC_TOLERANCE,
    ) -> "SymmetryContext":
        if not isinstance(setting, SpaceGroupSetting):
            raise TypeError("setting must be SpaceGroupSetting")
        provenance: tuple[tuple[str, object], ...] = (
            ("source", "space_group_setting"),
            ("setting_id", setting.setting_id),
            ("hall_symbol", setting.hall_symbol),
            ("catalog_number", setting.number),
        )
        return cls._build(
            setting.symmetry_operations,
            cell,
            setting_id=str(setting.setting_id),
            source_kind=SymmetrySourceKind.CATALOG_SETTING,
            diagnostics=(),
            provenance=provenance,
            metric_tolerance=metric_tolerance,
        )

    def operation_by_key(self, operation_key: str) -> AffineOperation:
        """Resolve a canonical operation key within this context."""

        try:
            return self._operation_lookup[operation_key]
        except KeyError as exc:
            raise KeyError(operation_key) from exc


__all__ = [
    "DEFAULT_METRIC_TOLERANCE",
    "DirectBasisConvention",
    "SymmetryContext",
    "SymmetryContextInvariantError",
    "SymmetrySourceKind",
    "canonical_operation_key",
]
