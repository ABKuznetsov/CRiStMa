from pathlib import Path

import cristma
from cristma.chemistry import (
    ChemistryAnalyzer,
    Composition,
    GrammarOperation,
    InteractionLayer,
)


FIXTURE = Path("tests/fixtures/cif/cod_3000098_barium_borate.cif")


def test_real_cif_reaches_actionable_chemistry_without_coordinates_in_analyzer() -> None:
    structure = cristma.read(FIXTURE).structures[0]

    result = ChemistryAnalyzer().analyze(Composition.from_structure(structure))

    assert result.composition.as_dict() == {"B": 36.0, "Ba": 18.0, "O": 72.0}
    assert result.classification.primary_family == "inorganic.oxide"
    assert {
        (
            interaction.centre_elements,
            interaction.operation,
            interaction.layer,
        )
        for interaction in result.grammar.candidate_interactions
    } == {
        (("B",), GrammarOperation.CENTRE_LIGAND_SHELL, InteractionLayer.STRUCTURAL),
        (("Ba",), GrammarOperation.INTERSTITIAL_COORDINATION, InteractionLayer.INTERSTITIAL),
    }
