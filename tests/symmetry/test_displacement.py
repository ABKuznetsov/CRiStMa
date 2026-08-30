import numpy as np
import pytest

from cristma.core.cell import UnitCell
from cristma.core.values import MeasuredValue
from cristma.structure import DisplacementParameters, IndependentSite, SiteComponent
from cristma.symmetry import (
    SymmetryConsistencyError,
    expand_orbit,
    parse_xyz_operation,
    transform_displacement,
)


def number(value: float, uncertainty: float | None = None) -> MeasuredValue:
    return MeasuredValue(value, uncertainty, str(value))


def anisotropic_u(values: tuple[tuple[float, float, float], ...]) -> DisplacementParameters:
    return DisplacementParameters(
        kind="U_aniso",
        tensor=tuple(tuple(number(value, 0.01) for value in row) for row in values),
        reported_kind="U",
    )


def numeric_tensor(displacement: DisplacementParameters) -> np.ndarray:
    assert displacement.tensor is not None
    return np.array(
        [[float(value.value) for value in row] for row in displacement.tensor],
        dtype=float,
    )


def test_anisotropic_u_tensor_rotates_with_symmetry_operation() -> None:
    displacement = anisotropic_u(
        ((1.0, 0.2, 0.0), (0.2, 2.0, 0.0), (0.0, 0.0, 3.0))
    )
    operation = parse_xyz_operation("-y,x,z", operation_id="op:rotate")

    transformed = transform_displacement(displacement, operation)

    rotation = np.asarray(operation.rotation, dtype=float)
    expected = rotation @ numeric_tensor(displacement) @ rotation.T
    assert np.allclose(numeric_tensor(transformed), expected)
    assert transformed.reported_kind == "U"
    assert transformed.tensor[0][0].raw is None


def test_adp_uncertainty_is_retained_when_all_contributing_terms_are_known() -> None:
    values = ((1.0, 0.2, 0.0), (0.2, 2.0, 0.0), (0.0, 0.0, 3.0))
    tensor = tuple(
        tuple(
            number(value, 0.01 if row == column else None)
            for column, value in enumerate(values[row])
        )
        for row in range(3)
    )
    displacement = DisplacementParameters(kind="U_aniso", tensor=tensor)

    transformed = transform_displacement(
        displacement,
        parse_xyz_operation("x,y,z", operation_id="op:identity"),
    )

    assert transformed.tensor[0][0].uncertainty == pytest.approx(0.01)
    assert transformed.tensor[0][1].uncertainty is None


def test_special_position_rejects_inconsistent_equivalent_adp_images() -> None:
    site = IndependentSite(
        id="site:X1",
        label="X1",
        components=(SiteComponent("Si", number(1.0)),),
        fractional=(number(0), number(0), number(0)),
        displacement=anisotropic_u(
            ((1.0, 0.0, 0.0), (0.0, 2.0, 0.0), (0.0, 0.0, 3.0))
        ),
    )
    operations = (
        parse_xyz_operation("x,y,z", operation_id="op:1"),
        parse_xyz_operation("y,x,z", operation_id="op:2"),
    )

    with pytest.raises(SymmetryConsistencyError, match="anisotropic displacement"):
        expand_orbit(site, operations, cell=UnitCell.cubic(number(4.0)))
