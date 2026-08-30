import numpy as np
import pytest

from cristma.core.cell import UnitCell
from cristma.core.structure import (
    Crystal,
    DisplacementParameters,
    IndependentSite,
    SiteComponent,
)
from cristma.core.values import MeasuredValue


def number(value: float) -> MeasuredValue:
    return MeasuredValue(value=value, uncertainty=None, raw=str(value))


def test_cell_rejects_non_positive_edge():
    with pytest.raises(ValueError, match="cell edge"):
        UnitCell(number(0), number(4), number(5), number(90), number(90), number(90))


def test_cubic_cell_exposes_metric_and_volume():
    cell = UnitCell.cubic(number(5.43))

    assert cell.volume == pytest.approx(5.43**3)
    assert np.diag(cell.metric) == pytest.approx([5.43**2] * 3)


def test_crystal_keeps_asymmetric_sites_primary():
    site = IndependentSite(
        id="site:Si1",
        label="Si1",
        components=(SiteComponent("Si", number(1.0)),),
        fractional=(number(0), number(0), number(0)),
    )
    crystal = Crystal(name="silicon", cell=UnitCell.cubic(number(5.43)), sites=(site,))

    assert crystal.sites == (site,)
    assert crystal.expanded_sites is None


def test_site_rejects_overoccupied_components():
    with pytest.raises(ValueError, match="occupancy"):
        IndependentSite(
            id="site:mixed",
            label="M1",
            components=(
                SiteComponent("La", number(0.7)),
                SiteComponent("Zr", number(0.4)),
            ),
            fractional=(number(0), number(0), number(0)),
        )


def test_site_keeps_reported_displacement_and_oxidation():
    displacement = DisplacementParameters(kind="U_iso", isotropic=number(0.01))
    site = IndependentSite(
        id="site:Zr1",
        label="Zr1",
        components=(SiteComponent("Zr", number(1), oxidation_state=number(4)),),
        fractional=(number(0), number(0), number(0)),
        displacement=displacement,
    )

    assert site.components[0].oxidation_state.value == 4
    assert site.displacement is displacement
