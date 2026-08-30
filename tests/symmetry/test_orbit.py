from cristma.core.structure import IndependentSite, SiteComponent
from cristma.core.values import MeasuredValue
from cristma.symmetry.affine import parse_xyz_operation
from cristma.symmetry.orbit import SpaceGroupDefinition, expand_orbit


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value=value, uncertainty=None, raw=str(value))


def site_at(x: float, y: float, z: float) -> IndependentSite:
    return IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1.0)),),
        fractional=(number(x), number(y), number(z)),
    )


def test_orbit_deduplicates_special_position_and_keeps_operation_ids():
    operations = (
        parse_xyz_operation("x,y,z", operation_id="op:1"),
        parse_xyz_operation("-x,-y,-z", operation_id="op:2"),
    )

    expanded = expand_orbit(site_at(0, 0, 0), operations)

    assert len(expanded) == 1
    assert expanded[0].independent_site_id == "site:Si1"
    assert expanded[0].equivalent_operation_ids == ("op:1", "op:2")


def test_orbit_keeps_translation_and_distinct_general_positions():
    operations = (
        parse_xyz_operation("x,y,z", operation_id="op:1"),
        parse_xyz_operation("-x,-y,-z", operation_id="op:2"),
    )

    expanded = expand_orbit(site_at(0.1, 0.2, 0.3), operations)

    assert [item.fractional for item in expanded] == [
        (0.1, 0.2, 0.3),
        (0.9, 0.8, 0.7),
    ]
    assert expanded[1].cell_translation == (-1, -1, -1)


def test_space_group_records_reported_provenance():
    identity = parse_xyz_operation("x,y,z", operation_id="op:1")

    group = SpaceGroupDefinition(
        hm_symbol="P 1",
        operations=(identity,),
        provenance="reported",
    )

    assert group.operations == (identity,)
    assert group.provenance == "reported"
