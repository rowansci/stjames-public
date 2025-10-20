from pytest import fixture, raises

from stjames import Mode, Molecule
from stjames.workflows import IRCWorkflow


@fixture
def water() -> Molecule:
    return Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")


def test_raise(water: Molecule) -> None:
    with raises(ValueError):
        IRCWorkflow(initial_molecule=water, mode=Mode.RECKLESS)


def test_smoke(water: Molecule) -> None:
    ircwf = IRCWorkflow(initial_molecule=water, mode=Mode.RAPID)

    assert ircwf.mode == Mode.RAPID
    assert ircwf.solvent is None
    assert not ircwf.preopt
    assert ircwf.max_irc_steps == 10
    assert ircwf.level_of_theory == "gfn2_xtb"
