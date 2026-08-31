import numpy as np
import pytest

from cristma.io.vasp.document import PoscarDocument, VaspScale, VaspSnapshot
from cristma.structure import SourceReference


def snapshot(**changes: object) -> VaspSnapshot:
    values = {
        "name": "demo",
        "lattice": np.eye(3),
        "species": ("Si",),
        "fractional": np.array([[0.25, 0.5, 0.75]]),
        "frame_index": 0,
        "source": SourceReference("POSCAR", "vasp", "frame:0", 0, 20),
    }
    values.update(changes)
    return VaspSnapshot(**values)


def test_document_preserves_complete_source() -> None:
    source = "demo\n1\n1 0 0\n0 1 0\n0 0 1\nSi\n1\nDirect\n0 0 0\n"

    document = PoscarDocument(raw_source=source, source_name="POSCAR")

    assert document.render_preserved() == source


def test_snapshot_normalizes_species_and_owns_read_only_arrays() -> None:
    original = np.array([[0.25, 0.5, 0.75]])
    value = snapshot(fractional=original)
    original[0, 0] = 0.9

    assert value.species[0].element == "Si"
    assert value.fractional.tolist() == [[0.25, 0.5, 0.75]]
    with pytest.raises(ValueError):
        value.fractional[0, 0] = 0.0


def test_velocity_requires_explicit_mode_and_unit() -> None:
    velocities = np.array([[1.0, 2.0, 3.0]])

    with pytest.raises(ValueError, match="velocity mode and unit"):
        snapshot(velocities=velocities)
    with pytest.raises(ValueError, match="velocity mode and unit"):
        snapshot(velocity_mode="cartesian", velocity_unit="angstrom/fs")

    value = snapshot(
        velocities=velocities,
        velocity_mode="cartesian",
        velocity_unit="angstrom/fs",
    )
    assert value.velocity_mode == "cartesian"


def test_snapshot_validates_atom_aligned_shapes_and_force_unit() -> None:
    with pytest.raises(ValueError, match="species"):
        snapshot(species=("Si", "O"))
    with pytest.raises(ValueError, match="selective_dynamics"):
        snapshot(selective_dynamics=np.ones((1, 2), dtype=bool))
    with pytest.raises(ValueError, match="force unit"):
        snapshot(forces=np.ones((1, 3)))


def test_scale_is_immutable_and_requires_supported_domain() -> None:
    assert VaspScale((1.0,)).values == (1.0,)

    with pytest.raises(ValueError):
        VaspScale((1.0, -2.0, 3.0))
