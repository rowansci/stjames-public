from typing import Any

from pytest import fixture, mark, raises

from stjames import Atom, Method, Mode, Molecule, Settings, SolventSettings, Task
from stjames.workflows import MultiStageOpt


@fixture
def He() -> Molecule:
    return Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])


@mark.parametrize(
    "mode, level_of_theory",
    [
        (Mode.RECKLESS, "gfn2_xtb/cpcmx(water)//gfn_ff"),
        (Mode.RAPID, "r2scan_3c/cpcm(water)//gfn2_xtb"),
        (Mode.CAREFUL, "wb97x_3c/cpcm(water)//b97_3c/cpcm(water)//gfn2_xtb"),
        (Mode.METICULOUS, "wb97m_d3bj/def2-tzvppd/cpcm(water)//wb97x_3c/cpcm(water)//b97_3c/cpcm(water)//gfn2_xtb"),
    ],
)
def test_multistage_opt_basic(mode: Mode, level_of_theory: str, He: Molecule) -> None:
    mso = MultiStageOpt(initial_molecule=He, mode=mode, solvent="water")

    assert str(mso) == f"<MultiStageOpt {mode}>"
    assert mso.level_of_theory == level_of_theory
    assert mso.optimization_settings
    assert mso.singlepoint_settings
    assert mso.singlepoint_settings.method
    assert mso.singlepoint_settings.solvent_settings
    assert mso.singlepoint_settings.solvent_settings.solvent == "water"
    assert mso.solvent == "water"
    assert not mso.constraints
    assert not mso.xtb_preopt
    assert not mso.transition_state


def test_raises(He: Molecule) -> None:
    with raises(ValueError):
        MultiStageOpt(initial_molecule=He)

    with raises(ValueError):
        # mypy is correct, but needs to be silenced to test pydantic error
        MultiStageOpt(mode="reckless", solvent="acetonitrile")  # type: ignore [call-arg]

    singlepoint_settings = Settings(method=Method.B973C)
    with raises(ValueError):
        MultiStageOpt(initial_molecule=He, mode=Mode.RAPID, singlepoint_settings=singlepoint_settings)

    with raises(NotImplementedError):
        MultiStageOpt(initial_molecule=He, mode="junk", solvent="acetonitrile")


def test_reckless(He: Molecule) -> None:
    mso = MultiStageOpt(initial_molecule=He, mode="reckless", solvent="acetonitrile")

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

    settings = {
        "method": Method.GFN_FF,
        "basis_set": None,
        "tasks": [Task.OPTIMIZE],
        "corrections": [],
        "mode": Mode.AUTO,
        "solvent_settings": None,
    }
    for key, value in settings.items():
        assert getattr(mso_opt0, key) == value

    assert not mso_opt0.opt_settings.transition_state


@mark.smoke
def test_rapid(He: Molecule) -> None:
    mso = MultiStageOpt(initial_molecule=He, mode="rapid", xtb_preopt=True, solvent="hexane")

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

    settings0 = {
        "method": Method.GFN0_XTB,
        "basis_set": None,
        "tasks": [Task.OPTIMIZE],
        "corrections": [],
        "mode": Mode.RAPID,
        "solvent_settings": None,
    }
    for key, value in settings0.items():
        assert getattr(mso_opt0, key) == value

    settings1 = {
        "method": Method.GFN2_XTB,
        "basis_set": None,
        "tasks": [Task.OPTIMIZE, Task.FREQUENCIES],
        "corrections": [],
        "mode": Mode.AUTO,
        "solvent_settings": None,
    }
    for key, value in settings1.items():
        assert getattr(mso_opt1, key) == value

    assert not mso_opt1.opt_settings.transition_state


def test_careful(He: Molecule) -> None:
    mso = MultiStageOpt(initial_molecule=He, mode="careful", solvent="acetone", transition_state=True)

    assert mso.optimization_settings
    assert len(mso.optimization_settings) == 2
    assert mso.singlepoint_settings
    assert mso.singlepoint_settings.method == Method.WB97X3C
    assert mso.singlepoint_settings.solvent_settings
    assert mso.singlepoint_settings.solvent_settings.solvent == "acetone"
    assert mso.singlepoint_settings.solvent_settings.model == "cpcm"
    assert mso.solvent == "acetone"
    assert not mso.constraints
    assert not mso.xtb_preopt
    assert mso.transition_state

    mso_opt0, mso_opt1 = mso.optimization_settings

    settings0: dict[str, Any] = {
        "method": Method.GFN2_XTB,
        "basis_set": None,
        "tasks": [Task.OPTIMIZE],
        "corrections": [],
        "mode": Mode.RAPID,
        "solvent_settings": None,
    }

    for key, value in settings0.items():
        assert getattr(mso_opt0, key) == value

    assert mso_opt0.opt_settings.transition_state

    settings1: dict[str, Any] = {
        "method": Method.B973C,
        "tasks": [Task.FREQUENCIES, Task.OPTIMIZE],
        "corrections": [],
        "mode": Mode.AUTO,
        "solvent_settings": SolventSettings(solvent="acetone", model="cpcm"),
    }
    for key, value in settings1.items():
        assert getattr(mso_opt1, key) == value

    assert mso_opt1.basis_set
    assert mso_opt1.basis_set.name == "def2-mTZVP"
    assert mso_opt1.opt_settings.transition_state


def test_meticulous(He: Molecule) -> None:
    mso = MultiStageOpt(initial_molecule=He, mode="meticulous", xtb_preopt=False)

    assert str(mso) == f"<MultiStageOpt {Mode.METICULOUS}>"
    assert mso.optimization_settings
    assert len(mso.optimization_settings) == 2
    assert mso.singlepoint_settings
    assert mso.singlepoint_settings.method == Method.WB97MD3BJ
    assert mso.singlepoint_settings.basis_set
    assert mso.singlepoint_settings.basis_set.name == "def2-TZVPPD"
    assert mso.singlepoint_settings.solvent_settings is None
    assert mso.solvent is None
    assert not mso.constraints
    assert not mso.xtb_preopt
    assert not mso.transition_state

    mso_opt0, mso_opt1 = mso.optimization_settings

    settings0 = {
        "method": Method.B973C,
        "tasks": [Task.OPTIMIZE],
        "corrections": [],
        "mode": Mode.AUTO,
        "solvent_settings": None,
    }

    for key, value in settings0.items():
        assert getattr(mso_opt0, key) == value

    assert mso_opt0.basis_set
    assert mso_opt0.basis_set.name == "def2-mTZVP"
    assert not mso_opt0.opt_settings.transition_state

    settings1 = {
        "method": Method.WB97X3C,
        "tasks": [Task.OPTIMIZE, Task.FREQUENCIES],
        "corrections": [],
        "mode": Mode.AUTO,
        "solvent_settings": None,
    }
    for key, value in settings1.items():
        assert getattr(mso_opt1, key) == value

    assert mso_opt1.basis_set
    assert mso_opt1.basis_set.name == "vDZP"
    assert not mso_opt1.opt_settings.transition_state


def test_manual(He: Molecule) -> None:
    """
    Builds a manual MultiStageOpt
    """
    optimization_settings = [
        Settings(method=Method.GFN0_XTB, tasks=[Task.OPTIMIZE]),
        Settings(method=Method.B3LYP, basis_set="def2-SVP", tasks=[Task.OPTIMIZE]),
    ]
    singlepoint_settings = Settings(method=Method.PBE, basis_set="def2-TZVP", solvent_settings=SolventSettings(solvent="octane", model="cpcm"))
    mso = MultiStageOpt(initial_molecule=He, optimization_settings=optimization_settings, singlepoint_settings=singlepoint_settings)

    level_of_theory = "pbe/def2-tzvp/cpcm(octane)//b3lyp/def2-svp//gfn0_xtb"

    assert str(mso) == f"<MultiStageOpt {level_of_theory}>"
    assert mso.level_of_theory == level_of_theory
    assert mso.optimization_settings
    assert len(mso.optimization_settings) == 2
    assert mso.singlepoint_settings
    assert mso.singlepoint_settings.method == Method.PBE
    assert mso.singlepoint_settings.basis_set
    assert mso.singlepoint_settings.basis_set.name == "def2-TZVP"
    assert mso.singlepoint_settings.solvent_settings
    assert mso.singlepoint_settings.solvent_settings.solvent == "octane"
    assert not mso.solvent
    assert not mso.constraints
    assert not mso.xtb_preopt
    assert not mso.transition_state

    mso_opt0, mso_opt1 = mso.optimization_settings

    settings0 = {
        "method": Method.GFN0_XTB,
        "tasks": [Task.OPTIMIZE],
        "corrections": [],
        "mode": Mode.AUTO,
        "solvent_settings": None,
    }

    for key, value in settings0.items():
        assert getattr(mso_opt0, key) == value

    assert not mso_opt0.opt_settings.transition_state

    settings1 = {
        "method": Method.B3LYP,
        "tasks": [Task.OPTIMIZE],
        "corrections": [],
        "mode": Mode.AUTO,
        "solvent_settings": None,
    }
    for key, value in settings1.items():
        assert getattr(mso_opt1, key) == value

    assert mso_opt1.basis_set
    assert mso_opt1.basis_set.name == "def2-SVP"
    assert not mso_opt1.opt_settings.transition_state
