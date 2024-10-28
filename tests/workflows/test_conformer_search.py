from pydantic.error_wrappers import ValidationError
from pytest import fixture, raises

from stjames import Constraint, Mode, Molecule
from stjames.method import Method
from stjames.workflows.conformer_search import (
    ConformerGenMixin,
    ConformerSearchMixin,
    ConformerSearchWorkflow,
    ETKDGSettings,
    ScreeningSettings,
    iMTDSettings,
    iMTDSpeeds,
)


@fixture
def water() -> Molecule:
    return Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")


@fixture
def ethane() -> Molecule:
    return Molecule.from_xyz("""\
C  -0.75  0.00  0.00
C   0.75  0.00  0.00
H  -1.15  0.00 -1.00
H  -1.15 -1.00  0.50
H  -1.15  0.85  0.50
H   1.15 -0.85 -0.50
H   1.15  0.85 -0.50
H   1.15  0.00  1.10
""")


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


def test_etkdg() -> None:
    reckless = ETKDGSettings(mode=Mode.RECKLESS)
    rapid = ETKDGSettings()

    assert str(reckless) == "<ETKDGSettings RECKLESS>"
    assert str(rapid) == "<ETKDGSettings RAPID>"
    assert reckless.num_initial_confs == 200
    assert rapid.num_initial_confs == 300

    with raises(NotImplementedError):
        ETKDGSettings(mode=Mode.CAREFUL)

    with raises(NotImplementedError):
        ETKDGSettings(mode=Mode.METICULOUS)

    with raises(ValidationError, match="ETKDG does not support NCI"):
        ETKDGSettings(nci=True)

    with raises(ValidationError, match="ETKDG does not support constraints"):
        ETKDGSettings(constraints=[Constraint(constraint_type="bond", atoms=[1, 2])])


def test_imtdgc() -> None:
    reckless = iMTDSettings(mode=Mode.RECKLESS)
    rapid = iMTDSettings()
    careful = iMTDSettings(mode=Mode.CAREFUL)
    meticulous = iMTDSettings(mode=Mode.METICULOUS)

    assert str(reckless) == "<iMTDSettings RECKLESS>"
    assert str(rapid) == "<iMTDSettings RAPID>"
    assert str(careful) == "<iMTDSettings CAREFUL>"
    assert str(meticulous) == "<iMTDSettings METICULOUS>"

    assert reckless.speed == iMTDSpeeds.MEGAQUICK
    assert rapid.speed == iMTDSpeeds.SUPERQUICK
    assert careful.speed == iMTDSpeeds.QUICK
    assert meticulous.speed == iMTDSpeeds.NORMAL

    assert reckless.reopt
    assert rapid.reopt
    assert not careful.reopt
    assert not meticulous.reopt

    assert reckless.mtd_method == Method.GFN_FF
    assert rapid.mtd_method == Method.GFN_FF
    assert careful.mtd_method == Method.GFN_FF
    assert meticulous.mtd_method == Method.GFN2_XTB

    assert reckless.conf_opt_method == Method.GFN_FF
    assert rapid.conf_opt_method == Method.GFN0_XTB
    assert careful.conf_opt_method == Method.GFN2_XTB
    assert meticulous.conf_opt_method == Method.GFN2_XTB


def test_conformer_gen_mixin() -> None:
    reckless = ConformerGenMixin(conf_gen_mode=Mode.RECKLESS)
    rapid = ConformerGenMixin()
    careful = ConformerGenMixin(conf_gen_mode=Mode.CAREFUL)
    meticulous = ConformerGenMixin(conf_gen_mode=Mode.METICULOUS)

    assert reckless.conf_gen_settings == ETKDGSettings(mode=Mode.RECKLESS)
    assert rapid.conf_gen_settings == ETKDGSettings()
    assert careful.conf_gen_settings == iMTDSettings(mode=Mode.CAREFUL)
    assert meticulous.conf_gen_settings == iMTDSettings(mode=Mode.METICULOUS)


def test_screening_settings() -> None:
    settings = ScreeningSettings(energy_threshhold=10)

    assert settings.energy_threshhold == 10
    assert settings.rotational_constants_threshhold == 0.02
    assert settings.max_confs is None


def test_conformer_search_mixin() -> None:
    reckless = ConformerSearchMixin(conf_gen_mode=Mode.RECKLESS, mso_mode=Mode.RECKLESS, solvent="water")
    rapid = ConformerSearchMixin(frequencies=True)
    careful_reckless = ConformerSearchMixin(conf_gen_mode=Mode.CAREFUL, mso_mode=Mode.RECKLESS, constraints=[Constraint(constraint_type="bond", atoms=[1, 2])])
    meticulous_rapid = ConformerSearchMixin(conf_gen_mode=Mode.METICULOUS, mso_mode=Mode.RAPID, screening_settings=ScreeningSettings(max_energy=100))

    assert reckless.conf_gen_mode == Mode.RECKLESS
    assert reckless.mso_mode == Mode.RECKLESS
    assert rapid.conf_gen_mode == Mode.RAPID
    assert rapid.mso_mode == Mode.RAPID
    assert careful_reckless.conf_gen_mode == Mode.CAREFUL
    assert careful_reckless.mso_mode == Mode.RECKLESS
    assert meticulous_rapid.conf_gen_mode == Mode.METICULOUS
    assert meticulous_rapid.mso_mode == Mode.RAPID

    assert reckless.solvent == "water"
    assert rapid.solvent is None

    assert rapid.frequencies
    assert not careful_reckless.frequencies

    assert careful_reckless.constraints == [Constraint(constraint_type="bond", atoms=[1, 2])]
    assert meticulous_rapid.constraints == ()


def test_conformer_search_workflow(water: Molecule) -> None:
    reckless_rapid = ConformerSearchWorkflow(initial_molecule=water, conf_gen_mode=Mode.RECKLESS, mso_mode=Mode.RAPID)
    rapid = ConformerSearchWorkflow(initial_molecule=water)
    careful = ConformerSearchWorkflow(initial_molecule=water, conf_gen_mode=Mode.CAREFUL, mso_mode=Mode.CAREFUL)
    meticulous = ConformerSearchWorkflow(
        initial_molecule=water,
        conf_gen_mode=Mode.METICULOUS,
        mso_mode=Mode.METICULOUS,
        constraints=[Constraint(constraint_type="bond", atoms=[1, 2])],
    )

    with raises(ValidationError, match="ETKDG does not support constraints"):
        ConformerSearchWorkflow(
            initial_molecule=water,
            conf_gen_mode=Mode.RECKLESS,
            constraints=[Constraint(constraint_type="bond", atoms=[1, 2])],
        )

    with raises(ValidationError, match="ETKDG does not support NCI"):
        ConformerSearchWorkflow(
            initial_molecule=water,
            conf_gen_mode=Mode.RECKLESS,
            nci=True,
        )

    assert reckless_rapid.conf_gen_mode == Mode.RECKLESS
    assert rapid.conf_gen_mode == Mode.RAPID
    assert careful.conf_gen_mode == Mode.CAREFUL
    assert meticulous.conf_gen_mode == Mode.METICULOUS

    assert reckless_rapid.mso_mode == Mode.RAPID
    assert rapid.mso_mode == Mode.RAPID
    assert careful.mso_mode == Mode.CAREFUL
    assert meticulous.mso_mode == Mode.METICULOUS

    assert reckless_rapid.conf_gen_settings == ETKDGSettings(mode=Mode.RECKLESS)
    assert rapid.conf_gen_settings == ETKDGSettings()
    assert careful.conf_gen_settings == iMTDSettings(mode=Mode.CAREFUL)
    assert meticulous.conf_gen_settings == iMTDSettings(mode=Mode.METICULOUS, constraints=[Constraint(constraint_type="bond", atoms=[1, 2])])
