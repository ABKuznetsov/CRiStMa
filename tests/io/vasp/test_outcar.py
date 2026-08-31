import numpy as np

from cristma.io.vasp.outcar import parse_outcar


TWO_IONIC_STEP_OUTCAR = """ vasp.6.4.3
 VRHFIN =Na: s1
 VRHFIN =Cl: s2p5
 ions per type = 1 1
 direct lattice vectors                 reciprocal lattice vectors
   2.0 0.0 0.0     0.5 0.0 0.0
   0.0 2.0 0.0     0.0 0.5 0.0
   0.0 0.0 2.0     0.0 0.0 0.5
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
   0.0 0.0 0.0      1.0 0.0 0.0
   1.0 1.0 1.0     -1.0 0.0 0.0
 DAV:  1  irrelevant electronic iteration
 direct lattice vectors                 reciprocal lattice vectors
   3.0 0.0 0.0     0.333 0.0 0.0
   0.0 3.0 0.0     0.0 0.333 0.0
   0.0 0.0 3.0     0.0 0.0 0.333
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----------------------------------------------------------------------------------
   0.0 0.0 0.0      0.1 0.2 0.3
   1.5 1.5 1.5     -0.1 -0.2 -0.3
"""


def test_outcar_uses_explicit_species_cell_positions_and_forces() -> None:
    result = parse_outcar(TWO_IONIC_STEP_OUTCAR, "OUTCAR")

    assert result.ok
    assert len(result.structures) == 2
    final = result.structures.final
    assert [site.components[0].element for site in final.sites] == ["Na", "Cl"]
    assert np.allclose(final.cell.matrix, np.diag([3.0, 3.0, 3.0]), atol=1e-12)
    assert final.properties["force"].unit == "eV/angstrom"
    assert final.atomic_view().properties["force"].values.shape == (2, 3)


def test_electronic_iterations_do_not_create_structure_frames() -> None:
    one_frame = TWO_IONIC_STEP_OUTCAR.split(" direct lattice vectors", 2)[0] + """
 direct lattice vectors                 reciprocal lattice vectors
   2 0 0  0 0 0
   0 2 0  0 0 0
   0 0 2  0 0 0
 DAV: 1
 DAV: 2
 POSITION                                       TOTAL-FORCE (eV/Angst)
 -----
 0 0 0  0 0 0
 1 1 1  0 0 0
 DAV: 3
"""

    result = parse_outcar(one_frame, "OUTCAR")

    assert len(result.structures) == 1


def test_truncated_final_position_block_is_ignored_with_diagnostic() -> None:
    source = TWO_IONIC_STEP_OUTCAR + """ POSITION TOTAL-FORCE (eV/Angst)
 -----
 0 0 0  0 0 0
"""

    result = parse_outcar(source, "OUTCAR")

    assert len(result.structures) == 2
    assert any(item.code == "vasp.outcar.frame_incomplete" for item in result.diagnostics)


def test_inconsistent_ions_per_type_prevents_frame_mapping() -> None:
    source = TWO_IONIC_STEP_OUTCAR.replace("ions per type = 1 1", "ions per type = 3")

    result = parse_outcar(source, "OUTCAR")

    assert not result.ok
    assert len(result.structures) == 0
    assert any(item.code == "vasp.outcar.species_count_inconsistent" for item in result.diagnostics)


def test_missing_species_are_explicit_unknowns() -> None:
    source = "\n".join(
        line for line in TWO_IONIC_STEP_OUTCAR.splitlines() if "VRHFIN" not in line
    )

    result = parse_outcar(source, "OUTCAR")

    assert result.ok
    assert result.structures.final.sites[0].components[0].element is None
    assert any(item.code == "vasp.map.species_unresolved" for item in result.diagnostics)


def test_titel_species_evidence_is_used_when_vrhfin_is_absent() -> None:
    source = TWO_IONIC_STEP_OUTCAR.replace(" VRHFIN =Na: s1", " TITEL = PAW_PBE Na_pv 08Apr2002").replace(
        " VRHFIN =Cl: s2p5", " TITEL = PAW_PBE Cl 06Sep2000"
    )

    result = parse_outcar(source, "OUTCAR")

    assert result.ok
    assert [site.components[0].element for site in result.structures.final.sites] == ["Na", "Cl"]
