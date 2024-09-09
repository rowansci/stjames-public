import pytest

from stjames import Atom, Molecule, MoleculeReadError, PeriodicCell

valid_extxyz = """
5
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

incorrect_num_atoms = """
6
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

not_digit_num_atoms = """
v
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

many_num_atoms = """
6 9
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""
no_num_atoms = """
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

xyz_style = """
5
Comment
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

missing_lattice = """
5
Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

missing_properties = """
5
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0"
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

incorrect_properites = """
5
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3foo:1
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

incorrect_lattice_extra = """
5
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0 3.14" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

incorrect_lattice_equals = """
5
Lattice="6.0 0.0 =0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

incorrect_lattice_str = """
5
Lattice="6.0 0.0 0.0 hi 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

incorrect_lattice_extra_string = """
5
Lattice="6.0 0.0 0.0 0.0 sup 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""


incorrect_lattice_single_quote = """
5
Lattice="6.0 0.0 0.0 0.0 6.0 '0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

incorrect_lattice_double_quote = """
5
Lattice="6.0 0.0 0.0 0.0 "6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

incorrect_lattice_double_single_quote = """
5
Lattice="6.0 0.0 0.0 0.0 '6.0 0.0 0.0 '0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

incorrect_lattice_double_double_quote = """
5
Lattice="6.0 0.0 "0.0 0.0 6.0 0.0 0.0 "0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""


expected_cell = (
    (6.0, 0.0, 0.0),
    (0.0, 6.0, 0.0),
    (0.0, 0.0, 6.0),
)

expected_atoms = [
    Atom(atomic_number=6, position=(0.0, 0.0, 0.0)),  # C
    Atom(atomic_number=1, position=(0.0, 0.0, 1.0)),  # H
    Atom(atomic_number=1, position=(1.0, 0.0, 0.0)),  # H
    Atom(atomic_number=1, position=(0.0, 1.0, 0.0)),  # H
    Atom(atomic_number=1, position=(1.0, 1.0, 1.0)),  # H
]

expected_molecule = Molecule(
    charge=0,
    multiplicity=1,
    atoms=expected_atoms,
    cell=PeriodicCell(lattice_vectors=expected_cell),
)


def test_molecule_from_extxyz_valid() -> None:
    """
    Test case for valid extxyz string.
    """
    molecule = Molecule.from_extxyz(valid_extxyz)
    assert molecule == expected_molecule, f"Valid case failed: got {molecule}, expected {expected_molecule}"


@pytest.mark.parametrize(
    "invalid_extxyz",
    [
        incorrect_num_atoms,
        no_num_atoms,
        not_digit_num_atoms,
        many_num_atoms,
        xyz_style,
        missing_lattice,
        missing_properties,
        incorrect_properites,
        incorrect_lattice_extra,
        incorrect_lattice_equals,
        incorrect_lattice_str,
        incorrect_lattice_extra_string,
        incorrect_lattice_single_quote,
        incorrect_lattice_double_quote,
        incorrect_lattice_double_single_quote,
        incorrect_lattice_double_double_quote,
    ],
)
def test_molecule_from_extxyz_invalid(invalid_extxyz: str) -> None:
    """
    Test case for invalid extxyz strings, ensuring they raise MoleculeReadError.
    """
    with pytest.raises(MoleculeReadError):
        Molecule.from_extxyz(invalid_extxyz)
