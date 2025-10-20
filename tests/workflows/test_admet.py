from pytest import raises

from stjames import ADMETWorkflow, Molecule


def test_input_type() -> None:
    molecule = Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")
    smiles = "c1ccccc1"

    ADMETWorkflow(initial_molecule=molecule)

    ADMETWorkflow(initial_smiles=smiles)

    with raises(ValueError):
        ADMETWorkflow(initial_molecule=molecule, initial_smiles=smiles)

    with raises(ValueError):
        ADMETWorkflow()
