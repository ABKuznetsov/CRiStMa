from pathlib import Path

import cristma


FIXTURE = Path(__file__).parents[2] / "fixtures" / "shelx" / "zdk288.res"


def test_real_shelx_fixture_maps_complete_crystal() -> None:
    result = cristma.read(FIXTURE)

    assert result.ok
    assert not result.diagnostics
    assert len(result.structures) == 1
    crystal = result.structures[0]
    assert crystal.name == "zdk288"
    assert len(crystal.sites) == 29
    assert len(crystal.space_group.operations) == 4
    assert len(crystal.atomic_view().atoms) == 116
    assert all(
        site.displacement is not None and site.displacement.kind == "U_aniso"
        for site in crystal.sites
    )


def test_real_shelx_source_round_trips_byte_for_byte(tmp_path: Path) -> None:
    result = cristma.read(FIXTURE)
    target = tmp_path / "zdk288-copy.res"

    cristma.write(result.document, target, mode="preserve")

    assert target.read_bytes() == FIXTURE.read_bytes()
