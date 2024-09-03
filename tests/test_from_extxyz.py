from stjames import Atom, Molecule, MoleculeReadError
import pytest


valid_extxyz = '''
5
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

incorrect_num_atoms = '''
6
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

no_num_atoms = '''
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

xyz_style = '''
5
Comment
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

missing_lattice =  '''
5
Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

missing_properties =  '''
5
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0"
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

incorrect_properites =  '''
5
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3foo:1
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

incorrect_lattice_extra =  '''
5
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0 3.14" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

incorrect_lattice_equals =  '''
5
Lattice="6.0 0.0 =0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

incorrect_lattice_str =  '''
5
Lattice="6.0 0.0 0.0 hi 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

incorrect_lattice_extra_string =  '''
5
Lattice="6.0 0.0 0.0 0.0 sup 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''


incorrect_lattice_single_quote =  '''
5
Lattice="6.0 0.0 0.0 0.0 6.0 '0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

incorrect_lattice_double_quote=  '''
5
Lattice="6.0 0.0 0.0 0.0 "6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

incorrect_lattice_double_single_quote =  '''
5
Lattice="6.0 0.0 0.0 0.0 '6.0 0.0 0.0 '0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''

incorrect_lattice_double_double_quote =  '''
5
Lattice="6.0 0.0 "0.0 0.0 6.0 0.0 0.0 "0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
'''


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
    cell=expected_cell,
)


def test_molecule_from_extxyz() -> None:

    test_cases = [
        ("valid_extxyz", valid_extxyz, True, expected_molecule),
        ("incorrect_num_atoms", incorrect_num_atoms, False, None),
        ("no_num_atoms", no_num_atoms, True, None),
        ("xyz_style", xyz_style, False, None),
        ("missing_lattice", missing_lattice, False, None),
        ("missing_properties", missing_properties, False, None),
        ("incorrect_properties", incorrect_properites, False, None),
        ("incorrect_lattice_extra", incorrect_lattice_extra, False, None),
        ("incorrect_lattice_equals", incorrect_lattice_equals, False, None),
        ("incorrect_lattice_str", incorrect_lattice_str, False, None),
        ("incorrect_lattice_extra_string", incorrect_lattice_extra_string, False, None),
        ("incorrect_lattice_single_quote", incorrect_lattice_single_quote, False, None),
        ("incorrect_lattice_double_quote", incorrect_lattice_double_quote, False, None),
        ("incorrect_lattice_double_single_quote", incorrect_lattice_double_single_quote, False, None),
    ]

    for name, extxyz, should_pass, expected in test_cases:
        print(f"========= TESTING {name}")
        if should_pass:
            molecule = Molecule.from_extxyz(extxyz)
            assert molecule == expected_molecule, f"Test {name} failed: molecule mismatch. Got {molecule} expected {expected}"
        else:
            try:
                with pytest.raises(MoleculeReadError):
                    Molecule.from_extxyz(extxyz)
            except AssertionError:
                pytest.fail(f"Test {name} failed: MoleculeReadError was not raised")

