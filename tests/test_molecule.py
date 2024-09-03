from functools import partial

from pydantic import ValidationError
from pytest import raises

from stjames import Atom, Molecule


def test_molecule_pbc() -> None:
    # this is a funny joke, because it makes manganese
    glomar_explorer = partial(Molecule, charge=0, multiplicity=2, atoms=[Atom(atomic_number=25, position=[0, 0, 0])])

    mol_nopbc = glomar_explorer()
    assert mol_nopbc.cell is None

    valid_input = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    mol = glomar_explorer(cell=valid_input)
    assert mol.cell == valid_input

    with raises(ValidationError):
        # Missing value
        glomar_explorer(cell=((1, 0, 0), (0, 1, 0), (0, 0)))

    with raises(ValidationError):
        # Bad type
        glomar_explorer(cell=((1, 0, 0), (0, 1, 0), (0, 0, "invalid")))

    with raises(ValidationError):
        # Too long
        glomar_explorer(cell=((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)))
