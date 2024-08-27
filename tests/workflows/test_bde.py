from pytest import fixture, raises

from stjames import Mode, Molecule
from stjames.workflows import BDEWorkflow


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


def test_raises(water: Molecule) -> None:
    with raises(ValueError):
        BDEWorkflow(initial_molecule=water, mode=Mode.RAPID, atoms=[5])


def test_ethane(ethane: Molecule) -> None:
    all_Hs = BDEWorkflow(initial_molecule=ethane, mode=Mode.RAPID, atoms=[3, 4, 5, 6, 7, 8])
    duplicated = BDEWorkflow(initial_molecule=ethane, mode=Mode.RAPID, atoms=[3, 4, 5, 6, 7, 8, 3, 4, 5, 6, 7, 8])
    all_CH = BDEWorkflow(initial_molecule=ethane, mode=Mode.RAPID, all_CH=True)
    all_CX = BDEWorkflow(initial_molecule=ethane, mode=Mode.RAPID, all_CX=True)
    ch3_frag_and_all_CH = BDEWorkflow(initial_molecule=ethane, mode=Mode.RAPID, fragments=[(2, 6, 7, 8)], all_CH=True)

    assert all_Hs.fragments == all_CH.fragments
    assert all_Hs.fragments == duplicated.fragments
    assert all_CX.fragments == ()
    assert ch3_frag_and_all_CH.fragments == ((2, 6, 7, 8), (3,), (4,), (5,), (6,), (7,), (8,))


def test_chloroethane(chloroethane: Molecule) -> None:
    all_Hs = BDEWorkflow(initial_molecule=chloroethane, mode=Mode.RAPID, atoms=[4, 5, 6, 7, 8])
    all_CH = BDEWorkflow(initial_molecule=chloroethane, mode=Mode.RAPID, all_CH=True)
    all_CX = BDEWorkflow(initial_molecule=chloroethane, mode=Mode.RAPID, all_CX=True)
    duplicated = BDEWorkflow(initial_molecule=chloroethane, mode=Mode.RAPID, all_CX=True, atoms=[3])
    ch3_frag_and_all_CH = BDEWorkflow(initial_molecule=chloroethane, mode=Mode.RAPID, fragments=[(2, 6, 7, 8)], all_CH=True)

    assert all_Hs.fragments == all_CH.fragments
    assert all_CX.fragments == ((3,),)
    assert duplicated.fragments == ((3,),)
    assert ch3_frag_and_all_CH.fragments == ((2, 6, 7, 8), (4,), (5,), (6,), (7,), (8,))
