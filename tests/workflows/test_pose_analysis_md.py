from pytest import fixture, raises

from stjames import Mode, Molecule
from stjames.pdb import PDB, read_pdb
from stjames.workflows import PoseAnalysisMolecularDynamicsWorkflow


@fixture
def water() -> Molecule:
    return Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")


@fixture
def gfp() -> PDB:
    """Green fluorescent protein."""
    return read_pdb("tests/data/1ema.pdb")


def test_raises(water: Molecule, gfp: str) -> None:
    PoseAnalysisMolecularDynamicsWorkflow(protein=gfp, initial_smiles="O", mode=Mode.RAPID)

    PoseAnalysisMolecularDynamicsWorkflow(protein=gfp, initial_smiles="O", mode=Mode.RAPID, protein_restraint_cutoff=5.0, protein_prune_cutoff=7.0)

    with raises(ValueError):
        PoseAnalysisMolecularDynamicsWorkflow(protein=gfp, initial_smiles="O", mode=Mode.RAPID, protein_restraint_cutoff=5.0, protein_prune_cutoff=3.0)
