from pydantic import ValidationError
from pytest import fixture, raises

from stjames import Constraint, Molecule
from stjames.conformers import ETKDGSettings, iMTDSettings
from stjames.workflows.conformer_search import ConformerSearchWorkflow


@fixture
def water() -> Molecule:
    return Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")


@fixture
def chloroethane() -> Molecule:
    return Molecule.from_xyz("""\
C  -0.75  0.00  0.00
C   0.75  0.00  0.00
Cl -1.15  0.00 -1.00
H  -1.15 -1.00  0.50
H  -1.15  0.85  0.50
H   1.15 -0.85 -0.50
H   1.15  0.85 -0.50
H   1.15  0.00  1.10
""")


def test_conformer_search_workflow(water: Molecule) -> None:
    ConformerSearchWorkflow(initial_molecule=water, conf_gen_settings=ETKDGSettings())
    ConformerSearchWorkflow(
        initial_molecule=water,
        conf_gen_settings=iMTDSettings(),
        constraints=[Constraint(constraint_type="bond", atoms=[1, 2])],
    )

    constraints = [Constraint(constraint_type="bond", atoms=[1, 2])]
    with raises(ValidationError, match="ETKDG does not support constraints"):
        ConformerSearchWorkflow(
            initial_molecule=water,
            conf_gen_settings=ETKDGSettings(constraints=constraints),
            constraints=constraints,
        )

    with raises(ValidationError, match="ETKDG does not support NCI"):
        ConformerSearchWorkflow(
            initial_molecule=water,
            conf_gen_settings=ETKDGSettings(nci=True),
            nci=True,
        )


def test_ts_constraints(water: Molecule) -> None:
    """Test that for transition states constraints are on for conformer generation and off for optimization."""
    constraints = [Constraint(constraint_type="bond", atoms=[1, 2])]
    imtd_settings = iMTDSettings(constraints=[Constraint(constraint_type="bond", atoms=[1, 2])])

    cr = ConformerSearchWorkflow(initial_molecule=water, constraints=constraints, transition_state=True, conf_gen_settings=imtd_settings)
    cr_msos = cr.multistage_opt_settings

    assert cr.conf_gen_settings is not None
    assert cr.conf_gen_settings.constraints == constraints
    assert len(cr_msos.optimization_settings) == 1

    assert not cr_msos.constraints
    assert cr_msos.optimization_settings[0].opt_settings.constraints == []


def test_no_conformer_gen(water: Molecule, chloroethane: Molecule) -> None:
    ConformerSearchWorkflow(initial_conformers=[chloroethane, chloroethane], conf_gen_settings=None)

    with raises(ValueError, match="Need initial_conformers to be set without a conformer-generation method"):
        ConformerSearchWorkflow(initial_molecule=water, conf_gen_settings=None)

    with raises(ValueError, match="Not all molecules in initial_conformers have the same atomic formula"):
        ConformerSearchWorkflow(initial_conformers=[chloroethane, water], conf_gen_settings=None)
