from pytest import fixture

from stjames import Mode, Molecule, Solvent
from stjames.workflows import RedoxPotentialWorkflow


@fixture
def water() -> Molecule:
    return Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")


def test_only_mecn(water: Molecule) -> None:
    wf = RedoxPotentialWorkflow(initial_molecule=water, mode=Mode.RAPID, solvent=Solvent.ACETONITRILE)
    assert wf.solvent == Solvent.ACETONITRILE

    wf = RedoxPotentialWorkflow(initial_molecule=water, mode=Mode.RAPID, solvent=None)
    assert wf.solvent == Solvent.ACETONITRILE

    wf = RedoxPotentialWorkflow(initial_molecule=water, mode=Mode.RAPID, solvent=Solvent.WATER)
    assert wf.solvent == Solvent.ACETONITRILE
