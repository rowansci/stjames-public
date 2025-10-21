from pytest import fixture, raises

from stjames import Molecule, pKaWorkflow


@fixture
def water() -> Molecule:
    return Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")


def test_chemprop_invalid_solvent(water: Molecule) -> None:
    """Test that ChemProp pKa prediction fails with an invalid solvent."""
    data = {
        "initial_smiles": "c1nc[nH]c1",
        "protonate_elements": [7],
        "deprotonate_elements": [7],
        "pka_range": (5, 15),
        "microscopic_pka_method": "chemprop_nevolianis2025",
        "solvent": "dimethylsulfoxide",
    }

    pKaWorkflow.model_validate(data)

    data2 = {
        "initial_molecule": water,
        "protonate_elements": [7],
        "deprotonate_elements": [7],
        "pka_range": (5, 15),
        "microscopic_pka_method": "aimnet2_wagen2024",
        "solvent": "dimethylsulfoxide",
    }

    with raises(ValueError):
        pKaWorkflow.model_validate(data2)

    data3 = {
        "initial_smiles": "c1nc[nH]c1",
        "protonate_elements": [7],
        "deprotonate_elements": [7],
        "pka_range": (5, 15),
        "microscopic_pka_method": "chemprop_nevolianis2025",
        "solvent": "invalid_solvent",
    }

    with raises(ValueError):
        pKaWorkflow.model_validate(data3)


def test_chemprop_extra_input(water: Molecule) -> None:
    """Test that ChemProp pKa prediction fails when initial_smiles and initial_molecule are given."""
    data = {
        "initial_molecule": water,
        "initial_smiles": "c1nc[nH]c1",
        "protonate_elements": [7],
        "deprotonate_elements": [7],
        "pka_range": (5, 15),
        "microscopic_pka_method": "chemprop_nevolianis2025",
    }

    with raises(ValueError):
        pKaWorkflow.model_validate(data)


def test_chemprop_missing_required_field() -> None:
    """Test that ChemProp pKa prediction fails when initial_smiles is missing."""
    data = {
        "protonate_elements": [7],
        "deprotonate_elements": [7],
        "pka_range": (5, 15),
        "microscopic_pka_method": "chemprop_nevolianis2025",
    }

    with raises(ValueError):
        pKaWorkflow.model_validate(data)


def test_chemprop_wrong_input_type() -> None:
    """Test that ChemProp pKa prediction fails with wrong input types."""
    data = {
        "initial_smiles": 12345,
        "protonate_elements": [7],
        "deprotonate_elements": [7],
        "pka_range": (5, 15),
        "microscopic_pka_method": "chemprop_nevolianis2025",
    }

    with raises(ValueError):
        pKaWorkflow.model_validate(data)


def test_chemprop_atoms_by_index_instead_of_element() -> None:
    """Test that using protonate_atoms/deprotonate_atoms with ChemProp will fail"""
    data = {
        "initial_smiles": "c1nc[nH]c1",
        "protonate_atoms": [1, 3],
        "deprotonate_atoms": [1, 3],
        "pka_range": (5, 15),
        "microscopic_pka_method": "chemprop_nevolianis2025",
    }

    with raises(ValueError):
        pKaWorkflow.model_validate(data)
