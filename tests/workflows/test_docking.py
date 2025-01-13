import atomium  # type: ignore [import-untyped]
from pytest import fixture, raises

from stjames import Mode, Molecule
from stjames.workflows import DockingWorkflow, Score


@fixture
def water() -> Molecule:
    return Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")


@fixture
def gfp() -> str:
    """Green fluorescent protein."""
    with open("tests/data/1ema.pdb") as f:
        return f.read()


def test_raises(water: Molecule, gfp: str) -> None:
    with raises(ValueError):
        DockingWorkflow(initial_molecule=water, mode=Mode.RAPID, target=gfp, pocket=((0, 0, 0), (-1, -1, -1)))


def test_basic(water: Molecule, gfp: str) -> None:
    dwf = DockingWorkflow(initial_molecule=water, mode=Mode.RAPID, target=gfp, pocket=((0, 0, 0), (1, 1, 1)))

    assert dwf.mode == Mode.RAPID
    assert not dwf.scores

    assert dwf.pocket == ((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))

    assert isinstance(dwf.target, atomium.data.File)
    assert dwf.target.code == "1EMA"
    assert dwf.target.filetype == "pdb"


def test_docked(water: Molecule, gfp: str) -> None:
    dwf = DockingWorkflow(initial_molecule=water, mode=Mode.RAPID, target=gfp, pocket=((0, 0, 0), (10, 10, 10)))

    pose1 = water.copy()
    pose2 = water.translated((1, 2, 3))

    dwf.scores = [
        Score(pose=pose1, score=0.0),
        Score(pose=pose2, score=1.0),
    ]

    assert dwf.pocket == ((0.0, 0.0, 0.0), (10.0, 10.0, 10.0))
    assert dwf.scores == [Score(pose=pose1, score=0.0), Score(pose=pose2, score=1.0)]
