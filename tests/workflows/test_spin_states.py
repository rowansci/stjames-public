from pytest import fixture, mark, raises

from stjames import Atom, Method, Mode, Molecule, Settings, Task
from stjames.workflows import MultiStageOptWorkflow, SpinStatesWorkflow


@fixture
def Mn() -> Molecule:
    return Molecule(charge=0, multiplicity=2, atoms=[Atom(atomic_number=25, position=[0, 0, 0])])


@fixture
def Fe() -> Molecule:
    return Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=26, position=[0, 0, 0])])


@mark.parametrize(
    "mode, level_of_theory",
    [
        (Mode.RECKLESS, "gfn2_xtb//gfn_ff"),
        (Mode.RAPID, "r2scan_3c//gfn2_xtb"),
        (Mode.CAREFUL, "wb97x_3c//b97_3c//gfn2_xtb"),
        (Mode.METICULOUS, "wb97m_d3bj/def2-tzvppd//wb97x_3c//b97_3c//gfn2_xtb"),
    ],
)
def test_spin_states_basic(mode: Mode, level_of_theory: str, Mn: Molecule) -> None:
    spin_states = SpinStatesWorkflow(
        initial_molecule=Mn,
        states=[2, 4, 6],
        mode=mode,
    )

    assert str(spin_states) == f"<SpinStatesWorkflow [2, 4, 6] {mode.name}>"
    assert spin_states.level_of_theory == level_of_theory
    assert spin_states.states == [2, 4, 6]

    msow = spin_states.multistage_opt_workflow
    assert msow
    assert msow.optimization_settings
    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method
    assert msow.singlepoint_settings.solvent_settings is None
    assert msow.solvent is None
    assert not msow.constraints
    assert msow.xtb_preopt is (mode in {Mode.CAREFUL, Mode.METICULOUS})
    assert not msow.transition_state


def test_raises(Mn: Molecule) -> None:
    with raises(ValueError):
        # mypy is correct, but needs to be silenced to test pydantic error
        SpinStatesWorkflow(initial_molecule=Mn, states=[1, 3, 5])  # type: ignore [call-arg]

    with raises(ValueError):
        # mypy is correct, but needs to be silenced to test pydantic error
        SpinStatesWorkflow(initial_molecule=Mn, mode=Mode.RECKLESS)  # type: ignore [call-arg]

    with raises(ValueError):
        SpinStatesWorkflow(
            initial_molecule=Mn,
            states=[1, 3, 4],
            mode=Mode.RECKLESS,
        )

    with raises(ValueError):
        SpinStatesWorkflow(
            initial_molecule=Mn,
            states=[2, 4, 5],
            mode=Mode.RAPID,
        )

    with raises(ValueError):
        SpinStatesWorkflow(
            initial_molecule=Mn,
            states=[2, 4, 6],
            mode=Mode.MANUAL,
        )


def test_auto(Mn: Molecule) -> None:
    spin_states = SpinStatesWorkflow(initial_molecule=Mn, mode=Mode.AUTO, states=[2, 4, 6], solvent="hexane")
    assert spin_states.mode == Mode.RAPID


def test_reckless(Mn: Molecule) -> None:
    spin_states = SpinStatesWorkflow(
        initial_molecule=Mn,
        states=[2, 4, 6],
        mode=Mode.RECKLESS,
        solvent="acetonitrile",
    )

    assert spin_states.states == [2, 4, 6]

    msow = spin_states.multistage_opt_workflow
    assert msow
    assert msow.optimization_settings
    assert len(msow.optimization_settings) == 1
    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method == Method.GFN2_XTB
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
def test_rapid(Mn: Molecule) -> None:
    spin_states = SpinStatesWorkflow(
        initial_molecule=Mn,
        states=[2, 4, 6],
        mode=Mode.RAPID,
        solvent="hexane",
        xtb_preopt=True,
    )

    assert spin_states.states == [2, 4, 6]

    msow = spin_states.multistage_opt_workflow
    assert msow
    assert msow.optimization_settings
    assert len(msow.optimization_settings) == 2
    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method == Method.R2SCAN3C
    assert msow.singlepoint_settings.solvent_settings
    assert msow.singlepoint_settings.solvent_settings.solvent == "hexane"
    assert msow.solvent == "hexane"
    assert not msow.constraints
    assert msow.xtb_preopt
    assert not msow.transition_state

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


def test_careful(Fe: Molecule) -> None:
    spin_states = SpinStatesWorkflow(
        initial_molecule=Fe,
        states=[1, 3, 5],
        mode=Mode.CAREFUL,
        transition_state=True,
    )

    assert spin_states.states == [1, 3, 5]

    msow = spin_states.multistage_opt_workflow
    assert msow
    assert msow.optimization_settings
    assert len(msow.optimization_settings) == 2
    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method == Method.WB97X3C
    assert not msow.singlepoint_settings.solvent_settings
    assert not msow.solvent
    assert not msow.constraints
    assert msow.xtb_preopt
    assert msow.transition_state

    msow_opt0, msow_opt1 = msow.optimization_settings

    assert msow_opt0.method == Method.GFN2_XTB
    assert msow_opt0.basis_set is None
    # Task.OPTIMIZE_TS is converted to Task.OPTIMIZE and settings.transition_state = True
    assert msow_opt0.tasks == [Task.OPTIMIZE]
    assert msow_opt0.corrections == []
    assert msow_opt0.mode == Mode.RAPID
    assert msow_opt0.solvent_settings is None
    assert msow_opt0.opt_settings.transition_state

    assert msow_opt1.method == Method.B973C
    assert msow_opt1.tasks == [Task.FREQUENCIES, Task.OPTIMIZE]
    assert msow_opt1.corrections == []
    assert msow_opt1.mode == Mode.AUTO
    assert msow_opt1.solvent_settings is None
    assert msow_opt1.basis_set
    assert msow_opt1.basis_set.name == "def2-mTZVP"
    assert msow_opt1.opt_settings.transition_state


def test_meticulous(Mn: Molecule) -> None:
    spin_states = SpinStatesWorkflow(
        initial_molecule=Mn,
        states=[2, 4, 6],
        mode=Mode.METICULOUS,
    )

    assert spin_states.states == [2, 4, 6]

    msow = spin_states.multistage_opt_workflow
    assert msow
    assert msow.optimization_settings
    assert len(msow.optimization_settings) == 3
    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method == Method.WB97MD3BJ
    assert msow.singlepoint_settings.basis_set
    assert msow.singlepoint_settings.basis_set.name == "def2-TZVPPD"
    assert msow.singlepoint_settings.solvent_settings is None
    assert msow.solvent is None
    assert not msow.constraints
    assert msow.xtb_preopt
    assert not msow.transition_state

    msow_opt0, msow_opt1, msow_opt2 = msow.optimization_settings

    assert msow_opt0.method == Method.GFN2_XTB
    assert msow_opt0.tasks == [Task.OPTIMIZE]
    assert msow_opt0.corrections == []
    assert msow_opt0.mode == Mode.RAPID
    assert msow_opt0.solvent_settings is None
    assert not msow_opt0.opt_settings.transition_state

    assert msow_opt1.method == Method.B973C
    assert msow_opt1.tasks == [Task.OPTIMIZE]
    assert msow_opt1.corrections == []
    assert msow_opt1.mode == Mode.AUTO
    assert msow_opt1.solvent_settings is None
    assert msow_opt1.basis_set
    assert msow_opt1.basis_set.name == "def2-mTZVP"
    assert not msow_opt1.opt_settings.transition_state

    assert msow_opt2.method == Method.WB97X3C
    assert msow_opt2.tasks == [Task.OPTIMIZE, Task.FREQUENCIES]
    assert msow_opt2.corrections == []
    assert msow_opt2.mode == Mode.AUTO
    assert msow_opt2.solvent_settings is None
    assert msow_opt2.basis_set
    assert msow_opt2.basis_set.name == "vDZP"
    assert not msow_opt2.opt_settings.transition_state


def test_manual(Fe: Molecule) -> None:
    """
    Builds a manual SpinStatesWorkflow
    """
    optimization_settings = [
        Settings(method=Method.GFN0_XTB, tasks=[Task.OPTIMIZE]),
        Settings(method=Method.B3LYP, basis_set="def2-SVP", tasks=[Task.OPTIMIZE]),
    ]
    singlepoint_settings = Settings(method=Method.PBE, basis_set="def2-TZVP")
    msow = MultiStageOptWorkflow(
        initial_molecule=Fe,
        mode=Mode.MANUAL,
        optimization_settings=optimization_settings,
        singlepoint_settings=singlepoint_settings,
    )
    spin_states = SpinStatesWorkflow(
        initial_molecule=Fe,
        mode=Mode.MANUAL,
        states=[1, 3, 5],
        multistage_opt_workflow=msow,
    )

    assert spin_states.states == [1, 3, 5]

    assert spin_states.multistage_opt_workflow
    msow = spin_states.multistage_opt_workflow
    assert msow
    assert msow.optimization_settings
    assert len(msow.optimization_settings) == 2
    assert msow.singlepoint_settings
    assert msow.singlepoint_settings.method == Method.PBE
    assert msow.singlepoint_settings.basis_set
    assert msow.singlepoint_settings.basis_set.name == "def2-TZVP"
    assert msow.singlepoint_settings.solvent_settings is None
    assert msow.solvent is None
    assert not msow.constraints
    assert not msow.xtb_preopt
    assert not msow.transition_state

    msow_opt0, msow_opt1 = msow.optimization_settings

    assert msow_opt0.method == Method.GFN0_XTB
    assert msow_opt0.tasks == [Task.OPTIMIZE]
    assert msow_opt0.corrections == []
    assert msow_opt0.mode == Mode.AUTO
    assert msow_opt0.solvent_settings is None
    assert not msow_opt0.opt_settings.transition_state

    assert msow_opt1.method == Method.B3LYP
    assert msow_opt1.tasks == [Task.OPTIMIZE]
    assert msow_opt1.corrections == []
    assert msow_opt1.mode == Mode.AUTO
    assert msow_opt1.solvent_settings is None

    assert msow_opt1.basis_set
    assert msow_opt1.basis_set.name == "def2-SVP"
    assert not msow_opt1.opt_settings.transition_state
