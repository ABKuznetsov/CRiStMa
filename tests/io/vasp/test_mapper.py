import numpy as np

from cristma.chemistry import UnknownSpecies
from cristma.io.vasp.document import VaspSnapshot
from cristma.io.vasp.mapper import map_vasp_snapshot
from cristma.structure import SourceReference


def test_snapshot_maps_cell_sites_and_typed_properties() -> None:
    source = SourceReference("CONTCAR", "vasp", "frame:2", 10, 90)
    snapshot = VaspSnapshot(
        name="relaxed",
        lattice=np.diag([2.0, 3.0, 4.0]),
        species=("Si", UnknownSpecies("vasp:type:2")),
        fractional=np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]),
        frame_index=2,
        source=source,
        selective_dynamics=np.array([[True, False, True], [False, False, False]]),
        velocities=np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]),
        velocity_mode="direct",
        velocity_unit="direct_lattice_vector/timestep",
        forces=np.array([[1.0, 2.0, 3.0], [-1.0, -2.0, -3.0]]),
        force_unit="eV/angstrom",
    )

    crystal = map_vasp_snapshot(snapshot)
    view = crystal.atomic_view()

    assert crystal.space_group.provenance == "unreported_identity"
    assert np.allclose(crystal.cell.matrix, snapshot.lattice, atol=1e-12)
    assert [site.components[0].element for site in crystal.sites] == ["Si", None]
    assert crystal.properties["selective_dynamics"].values.tolist() == [
        [True, False, True],
        [False, False, False],
    ]
    assert crystal.properties["velocity"].unit == "direct_lattice_vector/timestep"
    assert crystal.properties["velocity"].provenance.source_field == "velocity:direct"
    assert np.allclose(
        view.properties["force"].values,
        [[1.0, 2.0, 3.0], [-1.0, -2.0, -3.0]],
        atol=1e-12,
    )
    assert crystal.provenance.source is source


def test_cartesian_vectors_rotate_into_canonical_cell_frame() -> None:
    source_lattice = np.array(
        [
            [0.0, 2.0, 0.0],
            [-3.0, 0.0, 0.0],
            [0.0, 0.0, 4.0],
        ]
    )
    snapshot = VaspSnapshot(
        name="rotated",
        lattice=source_lattice,
        species=("Si",),
        fractional=np.array([[0.0, 0.0, 0.0]]),
        frame_index=0,
        source=SourceReference("OUTCAR", "vasp", "frame:0", 0, 10),
        velocities=np.array([[0.0, 1.0, 0.0]]),
        velocity_mode="cartesian",
        velocity_unit="angstrom/fs",
        forces=np.array([[0.0, 1.0, 0.0]]),
        force_unit="eV/angstrom",
    )

    crystal = map_vasp_snapshot(snapshot)

    assert np.allclose(crystal.cell.matrix, np.diag([2.0, 3.0, 4.0]), atol=1e-12)
    assert crystal.properties["velocity"].values.tolist() == [[1.0, 0.0, 0.0]]
    assert crystal.properties["force"].values.tolist() == [[1.0, 0.0, 0.0]]


def test_direct_velocities_are_not_cartesian_rotated() -> None:
    snapshot = VaspSnapshot(
        name="rotated",
        lattice=np.array([[0.0, 2.0, 0.0], [-3.0, 0.0, 0.0], [0.0, 0.0, 4.0]]),
        species=("Si",),
        fractional=np.array([[0.0, 0.0, 0.0]]),
        frame_index=0,
        source=SourceReference("CONTCAR", "vasp", "frame:0", 0, 10),
        velocities=np.array([[0.0, 1.0, 0.0]]),
        velocity_mode="direct",
        velocity_unit="direct_lattice_vector/timestep",
    )

    crystal = map_vasp_snapshot(snapshot)

    assert crystal.properties["velocity"].values.tolist() == [[0.0, 1.0, 0.0]]
