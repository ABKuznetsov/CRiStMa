import numpy as np
import pytest

from cristma.io.vasp.numeric import (
    fractional_from_cartesian,
    scaled_cartesian,
    scaled_lattice,
)
from cristma.io.vasp.document import VaspScale


def test_positive_scalar_scales_lattice_and_cartesian_positions() -> None:
    scale = VaspScale((2.0,))

    assert scaled_lattice(scale, np.eye(3)).tolist() == [
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.0, 0.0, 2.0],
    ]
    assert scaled_cartesian(scale, np.array([[1.0, 2.0, 3.0]]), np.eye(3)).tolist() == [
        [2.0, 4.0, 6.0]
    ]


def test_negative_scalar_reconstructs_requested_positive_volume() -> None:
    scale = VaspScale((-64.0,))
    lattice = scaled_lattice(scale, np.eye(3))

    assert np.linalg.det(lattice) == pytest.approx(64.0)
    assert scaled_cartesian(scale, np.array([[1.0, 0.0, 0.0]]), np.eye(3))[0, 0] == pytest.approx(4.0)


def test_three_scalars_scale_cartesian_components() -> None:
    scale = VaspScale((2.0, 3.0, 4.0))
    raw = np.array([[1.0, 1.0, 1.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]])

    assert scaled_lattice(scale, raw).tolist() == [
        [2.0, 3.0, 4.0],
        [0.0, 6.0, 0.0],
        [0.0, 0.0, 12.0],
    ]
    assert scaled_cartesian(scale, np.array([[1.0, 1.0, 1.0]]), raw).tolist() == [
        [2.0, 3.0, 4.0]
    ]


def test_cartesian_to_fractional_uses_rowwise_lattice_basis() -> None:
    lattice = np.diag([2.0, 4.0, 8.0])

    fractional = fractional_from_cartesian(lattice, np.array([[1.0, 2.0, 4.0]]))

    assert fractional.tolist() == [[0.5, 0.5, 0.5]]


@pytest.mark.parametrize(
    "values",
    [(0.0,), (-1.0, 2.0, 3.0), (1.0, 2.0), (float("nan"),)],
)
def test_invalid_scale_shape_or_domain_is_rejected(values: tuple[float, ...]) -> None:
    with pytest.raises(ValueError):
        VaspScale(values)


def test_negative_volume_rejects_singular_raw_lattice() -> None:
    with pytest.raises(ValueError, match="singular"):
        scaled_lattice(VaspScale((-10.0,)), np.zeros((3, 3)))
