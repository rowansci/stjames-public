from pytest import fixture, mark, raises

from stjames import Atom, Method, Mode, Molecule, Settings, SolventSettings, Task
from stjames.workflows import MultiStageOptWorkflow


@fixture
def He() -> Molecule:
    return Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])


@mark.parametrize(
    "mode, level_of_theory",
    [
        (Mode.RECKLESS, "gfn2_xtb/cpcmx(water)//gfn_ff/alpb(water)"),
        (Mode.RAPID, "r2scan_3c/cpcm(water)//gfn2_xtb/alpb(water)"),
        (Mode.CAREFUL, "wb97x_3c/cpcm(water)//r2scan_3c/cpcm(water)//gfn2_xtb"),
        (Mode.METICULOUS, "wb97m_d3bj/def2-tzvppd/cpcm(water)//wb97x_3c/cpcm(water)//r2scan_3c/cpcm(water)//gfn2_xtb"),
    ],
)
def test_multistage_opt_basic(mode: Mode, level_of_theory: str, He: Molecule) -> None:
    msow = MultiStageOptWorkflow(initial_molecule=He, mode=mode, solvent="water")

    assert str(msow) == f"<MultiStageOptWorkflow {mode.name}>"
    assert msow.level_of_theory == level_of_theory
    assert msow.optimization_settings
    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method
    assert msow.singlepoint_settings.solvent_settings
    assert msow.singlepoint_settings.solvent_settings.solvent == "water"
    assert msow.solvent == "water"
    assert not msow.constraints
    assert not msow.transition_state


def test_auto(He: Molecule) -> None:
    msow = MultiStageOptWorkflow(initial_molecule=He, mode=Mode.AUTO)
    assert msow.mode == Mode.RAPID


def test_raises(He: Molecule) -> None:
    with raises(ValueError):
        # mypy is correct, but needs to be silenced to test pydantic error
        MultiStageOptWorkflow(initial_molecule=He)  # type: ignore [call-arg]

    with raises(ValueError):
        # mypy is correct, but needs to be silenced to test pydantic error
        MultiStageOptWorkflow(mode=Mode.RAPID)  # type: ignore [call-arg]

    singlepoint_settings = Settings(method=Method.R2SCAN3C)
    with raises(ValueError):
        MultiStageOptWorkflow(initial_molecule=He, mode=Mode.RAPID, singlepoint_settings=singlepoint_settings)


def test_reckless(He: Molecule) -> None:
    msow = MultiStageOptWorkflow(initial_molecule=He, mode="reckless", solvent="acetonitrile")

    assert msow.optimization_settings
    assert len(msow.optimization_settings) == 1
    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method == Method.GFN2_XTB
    assert msow.singlepoint_settings.basis_set is None
    assert msow.singlepoint_settings.solvent_settings
    assert msow.singlepoint_settings.solvent_settings.solvent == "acetonitrile"
    assert msow.solvent == "acetonitrile"
    assert not msow.constraints
    assert not msow.xtb_preopt
    assert not msow.transition_state

    msow_opt0 = msow.optimization_settings[0]

    assert msow_opt0.method == Method.GFN_FF
    assert msow_opt0.basis_set is None
    assert msow_opt0.tasks == [Task.OPTIMIZE, Task.FREQUENCIES]
    assert msow_opt0.corrections == []
    assert msow_opt0.mode == Mode.AUTO
    assert msow_opt0.solvent_settings
    assert msow_opt0.solvent_settings.solvent == "acetonitrile"
    assert not msow_opt0.opt_settings.transition_state


@mark.smoke
def test_rapid(He: Molecule) -> None:
    msow = MultiStageOptWorkflow(initial_molecule=He, mode="rapid", xtb_preopt=True, solvent="hexane")

    assert msow.solvent == "hexane"
    assert not msow.constraints
    assert msow.xtb_preopt
    assert not msow.transition_state

    assert msow.optimization_settings
    assert len(msow.optimization_settings) == 2

    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method == Method.R2SCAN3C
    assert msow.singlepoint_settings.solvent_settings
    assert msow.singlepoint_settings.solvent_settings.solvent == "hexane"

    msow_opt0, msow_opt1 = msow.optimization_settings

    assert msow_opt0.method == Method.GFN0_XTB
    assert msow_opt0.basis_set is None
    assert msow_opt0.tasks == [Task.OPTIMIZE]
    assert msow_opt0.corrections == []
    assert msow_opt0.mode == Mode.RAPID
    assert msow_opt0.solvent_settings is None
    assert not msow_opt0.opt_settings.transition_state

    assert msow_opt1.method == Method.GFN2_XTB
    assert msow_opt1.basis_set is None
    assert msow_opt1.tasks == [Task.OPTIMIZE, Task.FREQUENCIES]
    assert msow_opt1.corrections == []
    assert msow_opt1.mode == Mode.AUTO
    assert msow_opt1.solvent_settings
    assert msow_opt1.solvent_settings.solvent == "hexane"
    assert not msow_opt1.opt_settings.transition_state


def test_careful(He: Molecule) -> None:
    msow = MultiStageOptWorkflow(initial_molecule=He, mode="careful", solvent="acetone", transition_state=True)

    assert msow.solvent == "acetone"
    assert not msow.constraints
    assert msow.xtb_preopt
    assert msow.transition_state

    assert msow.optimization_settings
    assert len(msow.optimization_settings) == 2

    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method == Method.WB97X3C
    assert msow.singlepoint_settings.solvent_settings
    assert msow.singlepoint_settings.solvent_settings.solvent == "acetone"
    assert msow.singlepoint_settings.solvent_settings.model == "cpcm"

    msow_opt0, msow_opt1 = msow.optimization_settings

    assert msow_opt0.method == Method.GFN2_XTB
    assert msow_opt0.basis_set is None
    assert msow_opt0.tasks == [Task.OPTIMIZE]
    assert msow_opt0.corrections == []
    assert msow_opt0.mode == Mode.RAPID
    assert msow_opt0.solvent_settings is None
    assert msow_opt0.opt_settings.transition_state

    assert msow_opt1.method == Method.R2SCAN3C
    assert msow_opt1.basis_set
    assert msow_opt1.basis_set.name == "def2-mTZVPP"
    assert msow_opt1.tasks == [Task.FREQUENCIES, Task.OPTIMIZE]
    assert msow_opt1.corrections == []
    assert msow_opt1.mode == Mode.AUTO
    assert msow_opt1.solvent_settings == SolventSettings(solvent="acetone", model="cpcm")
    assert msow_opt1.opt_settings.transition_state


def test_meticulous(He: Molecule) -> None:
    msow = MultiStageOptWorkflow(initial_molecule=He, mode="meticulous", xtb_preopt=False)

    assert msow.solvent is None
    assert not msow.constraints
    assert not msow.xtb_preopt
    assert not msow.transition_state
    assert str(msow) == "<MultiStageOptWorkflow METICULOUS>"

    assert msow.optimization_settings
    assert len(msow.optimization_settings) == 2

    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method == Method.WB97MD3BJ
    assert msow.singlepoint_settings.basis_set
    assert msow.singlepoint_settings.basis_set.name == "def2-TZVPPD"
    assert msow.singlepoint_settings.solvent_settings is None

    msow_opt0, msow_opt1 = msow.optimization_settings

    assert msow_opt0.method == Method.R2SCAN3C
    assert msow_opt0.basis_set
    assert msow_opt0.basis_set.name == "def2-mTZVPP"
    assert msow_opt0.tasks == [Task.OPTIMIZE]
    assert msow_opt0.corrections == []
    assert msow_opt0.mode == Mode.AUTO
    assert msow_opt0.solvent_settings is None
    assert not msow_opt0.opt_settings.transition_state

    assert msow_opt1.method == Method.WB97X3C
    assert msow_opt1.basis_set
    assert msow_opt1.basis_set.name == "vDZP"
    assert msow_opt1.tasks == [Task.OPTIMIZE, Task.FREQUENCIES]
    assert msow_opt1.corrections == []
    assert msow_opt1.mode == Mode.AUTO
    assert msow_opt1.solvent_settings is None
    assert not msow_opt1.opt_settings.transition_state


def test_manual(He: Molecule) -> None:
    """
    Builds a manual MultiStageOptWorkflow
    """
    optimization_settings = [
        Settings(method=Method.GFN0_XTB, tasks=[Task.OPTIMIZE]),
        Settings(method=Method.B3LYP, basis_set="def2-SVP", tasks=[Task.OPTIMIZE]),
    ]
    singlepoint_settings = Settings(method=Method.PBE, basis_set="def2-TZVP", solvent_settings=SolventSettings(solvent="octane", model="cpcm"))
    msow = MultiStageOptWorkflow(initial_molecule=He, mode=Mode.MANUAL, optimization_settings=optimization_settings, singlepoint_settings=singlepoint_settings)

    level_of_theory = "pbe/def2-tzvp/cpcm(octane)//b3lyp/def2-svp//gfn0_xtb"

    assert not msow.solvent
    assert not msow.constraints
    assert not msow.xtb_preopt
    assert not msow.transition_state
    assert str(msow) == f"<MultiStageOptWorkflow {level_of_theory}>"
    assert msow.level_of_theory == level_of_theory

    assert msow.optimization_settings
    assert len(msow.optimization_settings) == 2

    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method == Method.PBE
    assert msow.singlepoint_settings.basis_set
    assert msow.singlepoint_settings.basis_set.name == "def2-TZVP"
    assert msow.singlepoint_settings.solvent_settings
    assert msow.singlepoint_settings.solvent_settings.solvent == "octane"

    msow_opt0, msow_opt1 = msow.optimization_settings

    assert msow_opt0.method == Method.GFN0_XTB
    assert msow_opt0.basis_set is None
    assert msow_opt0.tasks == [Task.OPTIMIZE]
    assert msow_opt0.corrections == []
    assert msow_opt0.mode == Mode.AUTO
    assert msow_opt0.solvent_settings is None
    assert not msow_opt0.opt_settings.transition_state

    assert msow_opt1.method == Method.B3LYP
    assert msow_opt1.basis_set
    assert msow_opt1.basis_set.name == "def2-SVP"
    assert msow_opt1.tasks == [Task.OPTIMIZE]
    assert msow_opt1.corrections == []
    assert msow_opt1.mode == Mode.AUTO
    assert msow_opt1.solvent_settings is None
    assert not msow_opt1.opt_settings.transition_state
