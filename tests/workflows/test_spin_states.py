from pytest import fixture, mark, raises

from stjames import Atom, Method, Mode, Molecule, Settings, Task
from stjames.workflows import MultiStageOptInput, SpinStatesInput


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
    spin_states = SpinStatesInput(
        initial_molecule=Mn,
        states=[2, 4, 6],
        mode=mode,
    )

    assert str(spin_states) == f"<SpinStatesInput [2, 4, 6] {mode.name}>"
    assert spin_states.level_of_theory == level_of_theory
    assert spin_states.states == [2, 4, 6]

    mso = spin_states.multistage_opt_settings
    assert mso
    assert mso.optimization_settings
    assert mso.singlepoint_settings
    assert mso.singlepoint_settings.method
    assert mso.singlepoint_settings.solvent_settings is None
    assert mso.solvent is None
    assert not mso.constraints
    assert not mso.xtb_preopt
    assert not mso.transition_state


def test_raises(Mn: Molecule) -> None:
    with raises(ValueError):
        # mypy is correct, but needs to be silenced to test pydantic error
        SpinStatesInput(initial_molecule=Mn, states=[1, 3, 5])  # type: ignore [call-arg]

    with raises(ValueError):
        # mypy is correct, but needs to be silenced to test pydantic error
        SpinStatesInput(initial_molecule=Mn, mode=Mode.RECKLESS)  # type: ignore [call-arg]

    with raises(ValueError):
        SpinStatesInput(
            initial_molecule=Mn,
            states=[1, 3, 4],
            mode=Mode.RECKLESS,
        )

    with raises(ValueError):
        SpinStatesInput(
            initial_molecule=Mn,
            states=[2, 4, 5],
            mode=Mode.RAPID,
        )

    with raises(ValueError):
        SpinStatesInput(
            initial_molecule=Mn,
            states=[2, 4, 6],
            mode=Mode.MANUAL,
        )


def test_auto(Mn: Molecule) -> None:
    spin_states = SpinStatesInput(initial_molecule=Mn, mode=Mode.AUTO, states=[2, 4, 6], solvent="hexane")
    assert spin_states.mode == Mode.RAPID


def test_reckless(Mn: Molecule) -> None:
    spin_states = SpinStatesInput(
        initial_molecule=Mn,
        states=[2, 4, 6],
        mode=Mode.RECKLESS,
        solvent="acetonitrile",
    )

    assert spin_states.states == [2, 4, 6]

    mso = spin_states.multistage_opt_settings
    assert mso
    assert mso.optimization_settings
    assert len(mso.optimization_settings) == 1
    assert mso.singlepoint_settings
    assert mso.singlepoint_settings.method == Method.GFN2_XTB
    assert mso.singlepoint_settings.solvent_settings
    assert mso.singlepoint_settings.solvent_settings.solvent == "acetonitrile"
    assert mso.solvent == "acetonitrile"
    assert not mso.constraints
    assert not mso.xtb_preopt
    assert not mso.transition_state

    mso_opt0 = mso.optimization_settings[0]

    assert getattr(mso_opt0, "method") == Method.GFN_FF
    assert getattr(mso_opt0, "basis_set") is None
    assert getattr(mso_opt0, "tasks") == [Task.OPTIMIZE]
    assert getattr(mso_opt0, "corrections") == []
    assert getattr(mso_opt0, "mode") == Mode.AUTO
    assert getattr(mso_opt0, "solvent_settings", "Not None") is None
    assert not mso_opt0.opt_settings.transition_state


@mark.smoke
def test_rapid(Mn: Molecule) -> None:
    spin_states = SpinStatesInput(
        initial_molecule=Mn,
        states=[2, 4, 6],
        mode=Mode.RAPID,
        solvent="hexane",
        xtb_preopt=True,
    )

    assert spin_states.states == [2, 4, 6]

    mso = spin_states.multistage_opt_settings
    assert mso
    assert mso.optimization_settings
    assert len(mso.optimization_settings) == 2
    assert mso.singlepoint_settings
    assert mso.singlepoint_settings.method == Method.R2SCAN3C
    assert mso.singlepoint_settings.solvent_settings
    assert mso.singlepoint_settings.solvent_settings.solvent == "hexane"
    assert mso.solvent == "hexane"
    assert not mso.constraints
    assert mso.xtb_preopt
    assert not mso.transition_state

    mso_opt0, mso_opt1 = mso.optimization_settings

    assert mso_opt0.method == Method.GFN0_XTB
    assert mso_opt0.basis_set is None
    assert mso_opt0.tasks == [Task.OPTIMIZE]
    assert mso_opt0.corrections == []
    assert mso_opt0.mode == Mode.RAPID
    assert mso_opt0.solvent_settings is None
    assert not mso_opt0.opt_settings.transition_state

    assert mso_opt1.method == Method.GFN2_XTB
    assert mso_opt1.basis_set is None
    assert mso_opt1.tasks == [Task.OPTIMIZE, Task.FREQUENCIES]
    assert mso_opt1.corrections == []
    assert mso_opt1.mode == Mode.AUTO
    assert mso_opt1.solvent_settings is None
    assert not mso_opt1.opt_settings.transition_state


def test_careful(Fe: Molecule) -> None:
    spin_states = SpinStatesInput(
        initial_molecule=Fe,
        states=[1, 3, 5],
        mode=Mode.CAREFUL,
        transition_state=True,
    )

    assert spin_states.states == [1, 3, 5]

    mso = spin_states.multistage_opt_settings
    assert mso
    assert mso.optimization_settings
    assert len(mso.optimization_settings) == 2
    assert mso.singlepoint_settings
    assert mso.singlepoint_settings.method == Method.WB97X3C
    assert not mso.singlepoint_settings.solvent_settings
    assert not mso.solvent
    assert not mso.constraints
    assert not mso.xtb_preopt
    assert mso.transition_state

    mso_opt0, mso_opt1 = mso.optimization_settings

    assert mso_opt0.method == Method.GFN2_XTB
    assert mso_opt0.basis_set is None
    # Task.OPTIMIZE_TS is converted to Task.OPTIMIZE and settings.transition_state = True
    assert mso_opt0.tasks == [Task.OPTIMIZE]
    assert mso_opt0.corrections == []
    assert mso_opt0.mode == Mode.RAPID
    assert mso_opt0.solvent_settings is None
    assert mso_opt0.opt_settings.transition_state

    assert mso_opt1.method == Method.B973C
    assert mso_opt1.tasks == [Task.FREQUENCIES, Task.OPTIMIZE]
    assert mso_opt1.corrections == []
    assert mso_opt1.mode == Mode.AUTO
    assert mso_opt1.solvent_settings is None
    assert mso_opt1.basis_set
    assert mso_opt1.basis_set.name == "def2-mTZVP"
    assert mso_opt1.opt_settings.transition_state


def test_meticulous(Mn: Molecule) -> None:
    spin_states = SpinStatesInput(
        initial_molecule=Mn,
        states=[2, 4, 6],
        mode=Mode.METICULOUS,
    )

    assert spin_states.states == [2, 4, 6]

    mso = spin_states.multistage_opt_settings
    assert mso
    assert mso.optimization_settings
    assert len(mso.optimization_settings) == 3
    assert mso.singlepoint_settings
    assert mso.singlepoint_settings.method == Method.WB97MD3BJ
    assert mso.singlepoint_settings.basis_set
    assert mso.singlepoint_settings.basis_set.name == "def2-TZVPPD"
    assert mso.singlepoint_settings.solvent_settings is None
    assert mso.solvent is None
    assert not mso.constraints
    assert not mso.xtb_preopt
    assert not mso.transition_state

    mso_opt0, mso_opt1, mso_opt2 = mso.optimization_settings

    assert mso_opt0.method == Method.GFN2_XTB
    assert mso_opt0.tasks == [Task.OPTIMIZE]
    assert mso_opt0.corrections == []
    assert mso_opt0.mode == Mode.RAPID
    assert mso_opt0.solvent_settings is None
    assert not mso_opt0.opt_settings.transition_state

    assert mso_opt1.method == Method.B973C
    assert mso_opt1.tasks == [Task.OPTIMIZE]
    assert mso_opt1.corrections == []
    assert mso_opt1.mode == Mode.AUTO
    assert mso_opt1.solvent_settings is None
    assert mso_opt1.basis_set
    assert mso_opt1.basis_set.name == "def2-mTZVP"
    assert not mso_opt1.opt_settings.transition_state

    assert mso_opt2.method == Method.WB97X3C
    assert mso_opt2.tasks == [Task.OPTIMIZE, Task.FREQUENCIES]
    assert mso_opt2.corrections == []
    assert mso_opt2.mode == Mode.AUTO
    assert mso_opt2.solvent_settings is None
    assert mso_opt2.basis_set
    assert mso_opt2.basis_set.name == "vDZP"
    assert not mso_opt2.opt_settings.transition_state


def test_manual(Fe: Molecule) -> None:
    """
    Builds a manual SpinStatesInput
    """
    optimization_settings = [
        Settings(method=Method.GFN0_XTB, tasks=[Task.OPTIMIZE]),
        Settings(method=Method.B3LYP, basis_set="def2-SVP", tasks=[Task.OPTIMIZE]),
    ]
    singlepoint_settings = Settings(method=Method.PBE, basis_set="def2-TZVP")
    mso = MultiStageOptInput(
        initial_molecule=Fe,
        mode=Mode.MANUAL,
        optimization_settings=optimization_settings,
        singlepoint_settings=singlepoint_settings,
    )
    spin_states = SpinStatesInput(
        initial_molecule=Fe,
        mode=Mode.MANUAL,
        states=[1, 3, 5],
        multistage_opt_settings=mso,
    )

    assert spin_states.states == [1, 3, 5]

    assert spin_states.multistage_opt_settings
    mso = spin_states.multistage_opt_settings
    assert mso
    assert mso.optimization_settings
    assert len(mso.optimization_settings) == 2
    assert mso.singlepoint_settings
    assert mso.singlepoint_settings.method == Method.PBE
    assert mso.singlepoint_settings.basis_set
    assert mso.singlepoint_settings.basis_set.name == "def2-TZVP"
    assert mso.singlepoint_settings.solvent_settings is None
    assert mso.solvent is None
    assert not mso.constraints
    assert not mso.xtb_preopt
    assert not mso.transition_state

    mso_opt0, mso_opt1 = mso.optimization_settings

    assert mso_opt0.method == Method.GFN0_XTB
    assert mso_opt0.tasks == [Task.OPTIMIZE]
    assert mso_opt0.corrections == []
    assert mso_opt0.mode == Mode.AUTO
    assert mso_opt0.solvent_settings is None
    assert not mso_opt0.opt_settings.transition_state

    assert mso_opt1.method == Method.B3LYP
    assert mso_opt1.tasks == [Task.OPTIMIZE]
    assert mso_opt1.corrections == []
    assert mso_opt1.mode == Mode.AUTO
    assert mso_opt1.solvent_settings is None

    assert mso_opt1.basis_set
    assert mso_opt1.basis_set.name == "def2-SVP"
    assert not mso_opt1.opt_settings.transition_state
