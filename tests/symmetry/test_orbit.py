import numpy as np
import pytest

from cristma.core.cell import UnitCell
from cristma.core.structure import IndependentSite, SiteComponent
from cristma.core.values import MeasuredValue
from cristma.structure import (
    AtomicProperty,
    AtomicPropertyTable,
    CrystalStructure,
    PropertyProvenance,
    SymmetryImageProvenance,
)
from cristma.symmetry.affine import parse_xyz_operation
from cristma.symmetry.orbit import SpaceGroupDefinition, expand_orbit, expand_structure


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value=value, uncertainty=None, raw=str(value))


def site_at(x: float, y: float, z: float) -> IndependentSite:
    return IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1.0)),),
        fractional=(number(x), number(y), number(z)),
    )


def cubic_cell(edge: float = 4.0) -> UnitCell:
    return UnitCell.cubic(number(edge))


def test_orbit_deduplicates_special_position_and_keeps_operation_ids():
    operations = (
        parse_xyz_operation("x,y,z", operation_id="op:1"),
        parse_xyz_operation("-x,-y,-z", operation_id="op:2"),
    )

    expanded = expand_orbit(site_at(0, 0, 0), operations, cell=cubic_cell())

    assert len(expanded) == 1
    assert expanded[0].independent_site_id == "site:Si1"
    assert expanded[0].equivalent_images == (
        SymmetryImageProvenance("op:1", (0, 0, 0)),
        SymmetryImageProvenance("op:2", (0, 0, 0)),
    )
    assert expanded[0].representative_image in expanded[0].equivalent_images


def test_orbit_keeps_translation_and_distinct_general_positions():
    operations = (
        parse_xyz_operation("x,y,z", operation_id="op:1"),
        parse_xyz_operation("-x,-y,-z", operation_id="op:2"),
    )

    expanded = expand_orbit(site_at(0.1, 0.2, 0.3), operations, cell=cubic_cell())

    assert [item.fractional for item in expanded] == [
        (0.1, 0.2, 0.3),
        (0.9, 0.8, 0.7),
    ]
    assert expanded[1].representative_image.normalization_translation == (1, 1, 1)


def test_orbit_records_normalization_translation_with_explicit_sign() -> None:
    operation = parse_xyz_operation("x,y,z", operation_id="op:identity")

    atom = expand_orbit(
        site_at(1.12, -0.04, 0.35),
        (operation,),
        cell=cubic_cell(),
        structure_id="structure:test",
    )[0]

    assert atom.fractional == pytest.approx((0.12, 0.96, 0.35))
    assert atom.cartesian == pytest.approx((0.48, 3.84, 1.4))
    assert atom.representative_image == SymmetryImageProvenance(
        "op:identity", (-1, 1, 0)
    )


def test_orbit_keeps_mixed_occupation_on_one_geometric_position() -> None:
    site = IndependentSite(
        id="site:M1",
        label="M1",
        components=(
            SiteComponent("Ca", number(0.7)),
            SiteComponent("Sr", number(0.3)),
        ),
        fractional=(number(0), number(0), number(0)),
    )

    atoms = expand_orbit(
        site,
        (parse_xyz_operation("x,y,z", operation_id="op:1"),),
        cell=cubic_cell(),
    )

    assert len(atoms) == 1
    assert tuple(component.element for component in atoms[0].components) == ("Ca", "Sr")


def test_space_group_records_reported_provenance():
    identity = parse_xyz_operation("x,y,z", operation_id="op:1")

    group = SpaceGroupDefinition(
        hm_symbol="P 1",
        operations=(identity,),
        provenance="reported",
    )

    assert group.operations == (identity,)
    assert group.provenance == "reported"


def test_site_properties_reach_identity_atomic_view() -> None:
    operations = (parse_xyz_operation("x,y,z", operation_id="op:1"),)
    provenance = PropertyProvenance(source_name="OUTCAR", source_field="TOTAL-FORCE")
    crystal = CrystalStructure(
        "demo",
        cubic_cell(),
        (site_at(0.1, 0.2, 0.3),),
        id="structure:demo",
        space_group=SpaceGroupDefinition(operations, provenance="reported"),
        properties=AtomicPropertyTable(
            1,
            (
                AtomicProperty(
                    "force",
                    np.array([[1.0, 2.0, 3.0]]),
                    unit="eV/angstrom",
                    provenance=provenance,
                ),
            ),
        ),
    )

    expanded = expand_structure(crystal)

    assert expanded.properties["force"].values.tolist() == [
        [1.0, 2.0, 3.0],
    ]
    assert expanded.properties["force"].unit == "eV/angstrom"
    assert expanded.properties["force"].provenance is provenance


def test_nonidentity_property_expansion_requires_transform_semantics() -> None:
    operations = (
        parse_xyz_operation("x,y,z", operation_id="op:1"),
        parse_xyz_operation("-x,-y,-z", operation_id="op:2"),
    )
    crystal = CrystalStructure(
        "demo",
        cubic_cell(),
        (site_at(0.1, 0.2, 0.3),),
        space_group=SpaceGroupDefinition(operations, provenance="reported"),
        properties=AtomicPropertyTable(
            1,
            (AtomicProperty("force", np.array([[1.0, 2.0, 3.0]])),),
        ),
    )

    with pytest.raises(ValueError, match="transformation semantics"):
        expand_structure(crystal)
