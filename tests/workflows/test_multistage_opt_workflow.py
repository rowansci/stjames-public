from pytest import fixture, mark, raises

from stjames import Atom, Method, Mode, Molecule, Settings, Solvent, SolventSettings, Task
from stjames.workflows import MultiStageOptWorkflow, build_mso_settings, multi_stage_opt_settings_from_workflow


@fixture
def He() -> Molecule:
    return Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])


@mark.parametrize(
    "mode, level_of_theory",
    [
        (Mode.RECKLESS, "gfn2_xtb/cpcmx(water)//gfn_ff"),
        (Mode.RAPID, "r2scan_3c/cpcm(water)//gfn2_xtb"),
        (Mode.CAREFUL, "wb97x_3c/cpcm(water)//r2scan_3c//gfn2_xtb"),
        (Mode.METICULOUS, "wb97m_d3bj/def2-tzvppd/cpcm(water)//wb97x_3c//r2scan_3c//gfn2_xtb"),
    ],
)
def test_multistage_opt_basic(mode: Mode, level_of_theory: str, He: Molecule) -> None:
    msow = MultiStageOptWorkflow(initial_molecule=He, mode=mode, solvent="water", xtb_preopt=mode in {Mode.CAREFUL, Mode.METICULOUS})

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
        MultiStageOptWorkflow(mode=Mode.RAPID)  # type: ignore [call-arg]

    singlepoint_settings = Settings(method=Method.R2SCAN3C)
    with raises(ValueError):
        MultiStageOptWorkflow(initial_molecule=He, mode=Mode.RAPID, singlepoint_settings=singlepoint_settings)


def test_reckless(He: Molecule) -> None:
    msow = MultiStageOptWorkflow(initial_molecule=He, mode="reckless", solvent="acetonitrile", frequencies=True)

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
    assert msow_opt0.mode == Mode.RAPID
    assert msow_opt0.solvent_settings is None
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
    assert msow_opt1.tasks == [Task.OPTIMIZE]
    assert msow_opt1.corrections == []
    assert msow_opt1.mode == Mode.RAPID
    assert msow_opt1.solvent_settings is None
    assert not msow_opt1.opt_settings.transition_state


def test_careful(He: Molecule) -> None:
    msow = MultiStageOptWorkflow(initial_molecule=He, mode="careful", solvent="acetone", transition_state=True, xtb_preopt=True)

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
    assert msow_opt1.tasks == [Task.OPTIMIZE]
    assert msow_opt1.corrections == []
    assert msow_opt1.mode == Mode.RAPID
    assert msow_opt1.solvent_settings is None
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
    assert msow_opt0.mode == Mode.RAPID
    assert msow_opt0.solvent_settings is None
    assert not msow_opt0.opt_settings.transition_state

    assert msow_opt1.method == Method.WB97X3C
    assert msow_opt1.basis_set
    assert msow_opt1.basis_set.name == "vDZP"
    assert msow_opt1.tasks == [Task.OPTIMIZE]
    assert msow_opt1.corrections == []
    assert msow_opt1.mode == Mode.RAPID
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
    assert msow_opt0.mode == Mode.RAPID
    assert msow_opt0.solvent_settings is None
    assert not msow_opt0.opt_settings.transition_state

    assert msow_opt1.method == Method.B3LYP
    assert msow_opt1.basis_set
    assert msow_opt1.basis_set.name == "def2-SVP"
    assert msow_opt1.tasks == [Task.OPTIMIZE]
    assert msow_opt1.corrections == []
    assert msow_opt1.mode == Mode.RAPID
    assert msow_opt1.solvent_settings is None
    assert not msow_opt1.opt_settings.transition_state


def test_manual_from_factory(He: Molecule) -> None:
    """
    Builds a manual MultiStageOptWorkflow w/ new factory method
    """
    msos = build_mso_settings(
        sp_method=Method.PBE,
        sp_basis_set="def2-TZVP",
        opt_methods=[Method.GFN0_XTB, Method.B3LYP],
        opt_basis_sets=[None, "def2-SVP"],
        solvent=Solvent.OCTANE,
        use_solvent_for_opt=False,
    )

    level_of_theory = "pbe/def2-tzvp/cpcm(octane)//b3lyp/def2-svp//gfn0_xtb"

    assert not msos.constraints
    assert not msos.transition_state
    assert str(msos) == f"<MultiStageOptSettings {level_of_theory}>"
    assert msos.level_of_theory == level_of_theory

    assert msos.optimization_settings
    assert len(msos.optimization_settings) == 2

    assert msos.singlepoint_settings
    assert msos.singlepoint_settings.method == Method.PBE
    assert msos.singlepoint_settings.basis_set
    assert msos.singlepoint_settings.basis_set.name == "def2-TZVP"
    assert msos.singlepoint_settings.solvent_settings
    assert msos.singlepoint_settings.solvent_settings.solvent == "octane"

    msos_opt0, msos_opt1 = msos.optimization_settings

    assert msos_opt0.method == Method.GFN0_XTB
    assert msos_opt0.basis_set is None
    assert msos_opt0.tasks == [Task.OPTIMIZE]
    assert msos_opt0.corrections == []
    assert msos_opt0.mode == Mode.RAPID
    assert msos_opt0.solvent_settings is None
    assert not msos_opt0.opt_settings.transition_state

    assert msos_opt1.method == Method.B3LYP
    assert msos_opt1.basis_set
    assert msos_opt1.basis_set.name == "def2-SVP"
    assert msos_opt1.tasks == [Task.OPTIMIZE]
    assert msos_opt1.corrections == []
    assert msos_opt1.mode == Mode.RAPID
    assert msos_opt1.solvent_settings is None
    assert not msos_opt1.opt_settings.transition_state


def test_msos_from_msow(He: Molecule) -> None:
    msow = MultiStageOptWorkflow(initial_molecule=He, mode="careful", solvent="acetone", transition_state=True, xtb_preopt=True)

    msos = multi_stage_opt_settings_from_workflow(msow)

    assert msos.solvent == "acetone"
    assert not msos.constraints
    assert msos.xtb_preopt
    assert msos.transition_state

    assert msos.optimization_settings
    assert len(msos.optimization_settings) == 2

    assert msos.singlepoint_settings
    assert msos.singlepoint_settings.method == Method.WB97X3C
    assert msos.singlepoint_settings.solvent_settings
    assert msos.singlepoint_settings.solvent_settings.solvent == "acetone"
    assert msos.singlepoint_settings.solvent_settings.model == "cpcm"

    msos_opt0, msos_opt1 = msos.optimization_settings

    assert msos_opt0.method == Method.GFN2_XTB
    assert msos_opt0.basis_set is None
    assert msos_opt0.tasks == [Task.OPTIMIZE]
    assert msos_opt0.corrections == []
    assert msos_opt0.mode == Mode.RAPID
    assert msos_opt0.solvent_settings is None
    assert msos_opt0.opt_settings.transition_state

    assert msos_opt1.method == Method.R2SCAN3C
    assert msos_opt1.basis_set
    assert msos_opt1.basis_set.name == "def2-mTZVPP"
    assert msos_opt1.tasks == [Task.OPTIMIZE]
    assert msos_opt1.corrections == []
    assert msos_opt1.mode == Mode.RAPID
    assert msos_opt1.solvent_settings is None
    assert msos_opt1.opt_settings.transition_state
