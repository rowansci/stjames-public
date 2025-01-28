from pytest import fixture, raises

from stjames import Mode, Molecule
from stjames.pdb import PDB, read_pdb
from stjames.workflows import DockingWorkflow, Score


@fixture
def water() -> Molecule:
    return Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")


@fixture
def ammonia() -> Molecule:
    return Molecule.from_xyz("H 0 0 0\nN 0 0 1\nH 0 1 1\nH 1 0 1")


@fixture
def gfp() -> PDB:
    """Green fluorescent protein."""
    return read_pdb("tests/data/1ema.pdb")


def test_raises(water: Molecule, gfp: str) -> None:
    with raises(ValueError):
        DockingWorkflow(
            initial_molecule=water,
            mode=Mode.RAPID,
            smiles=["C", "N", "O"],
            target=gfp,
            pocket=((0, 0, 0), (-1, -1, -1)),
        )


def test_basic(water: Molecule, gfp: str) -> None:
    dwf = DockingWorkflow(
        initial_molecule=water,
        mode=Mode.RAPID,
        smiles=["C", "N", "O"],
        target=gfp,
        pocket=((0, 0, 0), (1, 1, 1)),
    )

    assert dwf.mode == Mode.RAPID
    assert not dwf.scores

    assert dwf.pocket == ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))

    assert isinstance(dwf.target, PDB)
    assert dwf.target.description.code == "1EMA"


def test_docked(water: Molecule, ammonia: Molecule, gfp: str) -> None:
    dwf = DockingWorkflow(
        initial_molecule=water,
        mode=Mode.RAPID,
        smiles=["N", "O"],
        target=gfp,
        pocket=((0, 0, 0), (10, 10, 10)),
    )

    dwf.scores = [
        Score(pose=water, score=0.0),
        Score(pose=ammonia, score=1.0),
    ]

    assert dwf.pocket == ((0.0, 0.0, 0.0), (10.0, 10.0, 10.0))
    assert dwf.scores == [Score(pose=water, score=0.0), Score(pose=ammonia, score=1.0)]
