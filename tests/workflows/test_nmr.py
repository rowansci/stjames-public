from pytest import fixture

from stjames import Mode, Molecule
from stjames.solvent import Solvent
from stjames.workflows import NMRSpectroscopyWorkflow


@fixture
def water() -> Molecule:
    return Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")


def test_basic(water: Molecule) -> None:
    nmr = NMRSpectroscopyWorkflow(
        initial_molecule=water,
        mode=Mode.RAPID,
    )

    assert nmr.mode == Mode.RAPID
    assert nmr.solvent == Solvent.CHLOROFORM

    nmr.per_conformer_isotropic_shieldings = [
        [0.0, 10.0, 20.0],
        [6.0, 12.0, 24.0],
    ]

    nmr.isotropic_shieldings = [3.0, 11.0, 22.0]

    assert len(nmr.isotropic_shieldings) == 3
    assert len(nmr.per_conformer_isotropic_shieldings) == 2
