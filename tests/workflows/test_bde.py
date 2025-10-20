from pytest import fixture, mark, raises

from stjames import Mode, Molecule
from stjames.workflows import BDE, BDEWorkflow


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


def test_auto(water: Molecule) -> None:
    wf = BDEWorkflow(initial_molecule=water, mode=Mode.AUTO, atoms=[1, 2])

    assert wf.mode == Mode.RAPID


def test_frequencies(water: Molecule) -> None:
    wf = BDEWorkflow(initial_molecule=water, mode=Mode.RAPID, atoms=[1, 2], frequencies=True)
    wf2 = BDEWorkflow(initial_molecule=water, mode=Mode.RAPID, atoms=[1, 2])
    wf3 = BDEWorkflow(initial_molecule=water, mode=Mode.RAPID, atoms=[1, 2], frequencies=False)

    assert wf.multistage_opt_settings.frequencies
    assert not wf2.multistage_opt_settings.frequencies
    assert not wf3.multistage_opt_settings.frequencies


def test_water(water: Molecule) -> None:
    all_Hs = BDEWorkflow(initial_molecule=water, mode=Mode.METICULOUS, atoms=[1, 2], optimize_fragments=True, xtb_preopt=True)
    duplicated = BDEWorkflow(initial_molecule=water, mode=Mode.METICULOUS, atoms=[1, 2, 1, 2])
    all_CH = BDEWorkflow(initial_molecule=water, mode=Mode.METICULOUS, all_CH=True, optimize_fragments=False, xtb_preopt=True)
    all_CX = BDEWorkflow(initial_molecule=water, mode=Mode.METICULOUS, all_CX=True)

    assert repr(all_Hs) == "<BDEWorkflow METICULOUS>"
    assert str(all_Hs) == "BDEWorkflow METICULOUS\n(1,)\n(2,)"

    assert all_Hs.level_of_theory == "wb97m_d3bj/def2-tzvppd//wb97x_3c//r2scan_3c//gfn2_xtb"

    assert all_Hs.fragment_indices == ((1,), (2,))
    assert duplicated.fragment_indices == ((1,), (2,))
    assert all_CH.fragment_indices == ()
    assert all_CX.fragment_indices == ()

    assert all_Hs.optimize_fragments
    assert duplicated.optimize_fragments
    assert not all_CH.optimize_fragments
    assert all_CX.optimize_fragments


def test_ethane(ethane: Molecule) -> None:
    all_Hs = BDEWorkflow(initial_molecule=ethane, mode=Mode.RAPID, atoms=[3, 4, 5, 6, 7, 8], optimize_fragments=True)
    duplicated = BDEWorkflow(initial_molecule=ethane, mode=Mode.RAPID, atoms=[3, 4, 5, 6, 7, 8, 3, 4, 5, 6, 7, 8])
    all_CH = BDEWorkflow(initial_molecule=ethane, mode=Mode.RAPID, all_CH=True, optimize_fragments=False)
    all_CX = BDEWorkflow(initial_molecule=ethane, mode=Mode.RAPID, all_CX=True)
    ch3_frag_and_all_CH = BDEWorkflow(initial_molecule=ethane, mode=Mode.RAPID, fragment_indices=[(2, 6, 7, 8)], all_CH=True)

    assert repr(all_Hs) == "<BDEWorkflow RAPID>"
    assert str(all_Hs) == "BDEWorkflow RAPID\n(3,)\n(4,)\n(5,)\n(6,)\n(7,)\n(8,)"

    assert all_Hs.fragment_indices == all_CH.fragment_indices
    assert all_Hs.fragment_indices == duplicated.fragment_indices
    assert all_CX.fragment_indices == ()
    assert ch3_frag_and_all_CH.fragment_indices == ((2, 6, 7, 8), (3,), (4,), (5,), (6,), (7,), (8,))

    assert all_Hs.optimize_fragments
    assert duplicated.optimize_fragments
    assert not all_CH.optimize_fragments
    assert all_CX.optimize_fragments
    assert ch3_frag_and_all_CH.optimize_fragments


def test_chloroethane(chloroethane: Molecule) -> None:
    all_Hs = BDEWorkflow(initial_molecule=chloroethane, mode=Mode.RAPID, atoms=[4, 5, 6, 7, 8])
    all_CH = BDEWorkflow(initial_molecule=chloroethane, mode=Mode.RAPID, all_CH=True)
    all_CX = BDEWorkflow(initial_molecule=chloroethane, mode=Mode.RAPID, all_CX=True)
    duplicated = BDEWorkflow(initial_molecule=chloroethane, mode=Mode.RAPID, all_CX=True, atoms=[3])
    ch3_frag_and_all_CH = BDEWorkflow(initial_molecule=chloroethane, mode=Mode.RAPID, fragment_indices=[(2, 6, 7, 8)], all_CH=True)

    assert repr(all_Hs) == "<BDEWorkflow RAPID>"
    assert str(all_Hs) == "BDEWorkflow RAPID\n(4,)\n(5,)\n(6,)\n(7,)\n(8,)"

    assert all_Hs.fragment_indices == all_CH.fragment_indices
    assert all_CX.fragment_indices == ((3,),)
    assert duplicated.fragment_indices == ((3,),)
    assert ch3_frag_and_all_CH.fragment_indices == ((2, 6, 7, 8), (4,), (5,), (6,), (7,), (8,))


@mark.parametrize(
    "mode, opt_frag",
    [
        (Mode.RECKLESS, False),
        (Mode.RAPID, True),
        (Mode.CAREFUL, True),
        (Mode.METICULOUS, True),
    ],
)
def test_mode_defaults(chloroethane: Molecule, mode: Mode, opt_frag: bool) -> None:
    wf = BDEWorkflow(initial_molecule=chloroethane, mode=mode, atoms=[4, 5, 6, 7, 8])

    assert wf.fragment_indices == ((4,), (5,), (6,), (7,), (8,))
    assert wf.optimize_fragments == opt_frag


def test_BDE_none() -> None:
    """Test that a BDE with None for an energy is fine."""

    BDE(fragment_idxs=(1, 2), energy=None, fragment_energies=(1, None), calculation_uuids=(["1", "2"], ["2", "3"]))
    BDE(fragment_idxs=(1, 2), energy=None, fragment_energies=(None, 2), calculation_uuids=(["1", "2"], ["2", "3"]))
    BDE(fragment_idxs=(1, 2), energy=None, fragment_energies=(1, 2), calculation_uuids=(["1", "2"], ["2", "3"]))


def test_new_methods(chloroethane: Molecule) -> None:
    g_xtb = BDEWorkflow(initial_molecule=chloroethane, mode="G_XTB", atoms=[4, 5, 6, 7, 8])
    g_xtb__gfn2_xtb = BDEWorkflow(initial_molecule=chloroethane, mode="G_XTB//GFN2_XTB", all_CH=True)
    BDEWorkflow(initial_molecule=chloroethane, mode="R2SCAN3C//GFN2_XTB", all_CX=True)
    BDEWorkflow(initial_molecule=chloroethane, mode="OMOL25_CONSERVING_S", all_CX=True, atoms=[3])

    assert str(g_xtb) == "BDEWorkflow G_XTB\n(4,)\n(5,)\n(6,)\n(7,)\n(8,)"
    assert repr(g_xtb__gfn2_xtb) == "<BDEWorkflow G_XTB//GFN2_XTB>"
