from functools import partial

from pydantic import ValidationError
from pytest import raises

from stjames import Atom, Molecule


def test_molecule_pbc() -> None:
    # this is a funny joke, because it makes manganese
    glomar_explorer = partial(Molecule, charge=0, multiplicity=2, atoms=[Atom(atomic_number=25, position=[0, 0, 0])])

    mol_nopbc = glomar_explorer()
    assert mol_nopbc.cell is None

    valid_input = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    mol = glomar_explorer(cell={"lattice_vectors": valid_input})
    assert mol.cell is not None
    assert mol.cell.lattice_vectors == valid_input

    mol2 = glomar_explorer(cell={"lattice_vectors": valid_input, "is_periodic": (True, False, True)})
    assert mol2.cell is not None
    assert mol2.cell.lattice_vectors == valid_input
    assert mol2.cell.is_periodic == (True, False, True)

    with raises(ValidationError):
        glomar_explorer(cell={"lattice_vectors": valid_input, "is_periodic": (False, False, False)})

    invalid_input = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0))
    with raises(ValidationError):
        glomar_explorer(cell={"lattice_vectors": invalid_input})

    invalid_input2 = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, "invalid"))
    with raises(ValidationError):
        glomar_explorer(cell={"lattice_vectors": invalid_input2})

    invalid_input3 = ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 1.0, 1.0))
    with raises(ValidationError):
        glomar_explorer(cell={"lattice_vectors": invalid_input3})


def test_xyz() -> None:
    mol = Molecule.from_xyz_lines(
        """3
H2O Molecule; smiles: O; nonsense: data: more nonsense
O        0.000000     0.000000     0.000000
H        0.000000    -0.757000     0.587000
H        0.000000     0.757000     0.587000
""".splitlines()
    )
    assert mol.smiles == "O"
    assert getattr(mol, "nonsense", None) is None

    mol = Molecule.from_xyz_lines(
        """3
name: H2O Molecule; smiles: O
O        0.000000     0.000000     0.000000
H        0.000000    -0.757000     0.587000
H        0.000000     0.757000     0.587000
""".splitlines()
    )

    assert getattr(mol, "name", None) is None

    mol = Molecule.from_xyz_lines(
        """3
cell: [(1, 2, 3), [4, 5, 6e-1], (7.0, 8.0, 9.0)]; smiles: O; is_periodic: [True, False, True]
O        0.000000     0.000000     0.000000
H        0.000000    -0.757000     0.587000
H        0.000000     0.757000     0.587000
""".splitlines()
    )

    assert mol.cell is not None
    assert mol.cell.lattice_vectors == ((1.0, 2.0, 3.0), (4.0, 5.0, 0.6), (7.0, 8.0, 9.0))
    assert mol.cell.is_periodic == (True, False, True)
    assert mol.smiles == "O"

    mol = Molecule.from_xyz_lines(
        """3
cell: [(a, b, c), [4, 5, 6e-1], (7.0, 8.0, 9.0)]; smiles: O; is_periodic: [True, False, True]
O        0.000000     0.000000     0.000000
H        0.000000    -0.757000     0.587000
H        0.000000     0.757000     0.587000
""".splitlines()
    )

    assert mol.cell is None
    assert mol.smiles == "O"
    assert getattr(mol.cell, "is_periodic", None) is None
