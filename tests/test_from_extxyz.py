from pytest import mark, raises

from stjames import Atom, Molecule, MoleculeReadError, PeriodicCell

# Valid EXTXYZ without forces (only positions)
valid_extxyz = """
5
Lattice="6.0 0.0 0.0 0.0 6.0 0.0 0.0 0.0 6.0" Properties=species:S:1:pos:R:3
C        0.0        0.0        0.0
H        0.0        0.0        1.0
H        1.0        0.0        0.0
H        0.0        1.0        0.0
H        1.0        1.0        1.0
"""

# Valid EXTXYZ with forces: forces are provided after positions.
# For each atom, the forces are given and the molecule's gradient should be -forces.
valid_extxyz_with_forces_and_energy = """
22
Properties=species:S:1:pos:R:3:forces:R:3:MACE_forces:R:3 smiles="[H:13][C:5]([H:14])([H:15])[C:7]([H:18])([H:19])[C:6]([H:16])([H:17])[S:2][C:4]([H:11])([H:12])[S:1][C:3]([H:8])([H:9])[H:10].[H:21][O:20][H:22]" total_charge=0 energy=-29139.55729148458 config_type="DES370K Dimers" MACE_energy=-29139.53445598132 pbc="F F F"
S        2.75637007       2.70758009       4.82925987       1.55959951      -0.22801913      -0.54136018       1.57411973      -0.23083920      -0.56575688
S        2.84554004       0.22606000       3.09241009      -0.21334177      -1.76227773       0.33161567      -0.21053821      -1.76733183       0.31476973
C        1.01258004       2.49340009       4.93010997      -1.31854381      -0.37388875      -0.24927582      -1.32686027      -0.41468579      -0.24419847
C        3.17725992       1.94638002       3.18563008      -0.09910375       1.52775439      -0.55371287      -0.09844947       1.52902410      -0.49063227
C        0.00000000       0.00000000       0.00000000       0.54047842      -1.02932160      -0.37462899       0.55420524      -1.04859973      -0.38165141
C        2.22169995      -0.06399999       1.41622007      -0.74959058       0.13849016      -0.66947629      -0.74189889       0.14272211      -0.65926828
C        0.72196996       0.27478000       1.31911993      -0.88825228       1.14448169       0.10462363      -0.88035930       1.15971257       0.11942591
H        0.73436999       1.43289995       4.80946970       0.11871475       0.45471041       0.08255662       0.12983545       0.47803891       0.09119142
H        0.46778000       3.08784985       4.17749977       0.24911941      -0.20360289       0.32966949       0.24533754      -0.20850576       0.33744960
H        0.70124000       2.82674027       5.93438053       0.20534666      -0.13940533      -0.46775994       0.20706729      -0.12047003      -0.47116810
H        2.59819007       2.49752998       2.41988992       0.29845053      -0.31965051       0.34195845       0.29323885      -0.33151445       0.32609726
H        4.25116014       2.06114006       2.94268989      -0.62592601       0.53239934       0.63500463      -0.63988801       0.53433737       0.61779540
H        0.01781000      -1.10166001      -0.05143000      -0.31652062       0.33261293      -0.62734778      -0.32052002       0.33247610      -0.62978105
H        0.50296998       0.40217999      -0.89547998      -0.14484048      -0.10780390       0.42454441      -0.15022920      -0.10846825       0.42265847
H       -1.05070996       0.33554000      -0.00612000       0.55084403      -0.10031390       0.18450139       0.55970312      -0.09816881       0.19243150
H        2.77840995       0.53621000       0.67106002      -0.14551747      -0.35971650       0.43068238      -0.14826788      -0.34813759       0.42951233
H        2.39210987      -1.13260996       1.18280005       0.04825265       0.56114332       0.26057253       0.04712810       0.55862479       0.26730085
H        0.59988999       1.35992002       1.50082004       0.18014249      -0.45468505       0.36480505       0.17898520      -0.45518615       0.36105601
H        0.10381000      -0.26278001       2.06366992       0.77321355       0.35561432       0.04327458       0.76235138       0.36555290       0.03043141
O        5.42839002       0.27989000       7.42625999      -0.71867266       0.97644353       0.67783188      -0.72326306       0.98088507       0.67544292
H        4.56538010       0.65290999       7.13586997       1.00844588      -0.37843304       0.30550308       0.99937753      -0.37383106       0.29632339
H        5.64198017       0.76508999       8.25524997      -0.31126251      -0.56856677      -1.03720163      -0.31107513      -0.57563526      -1.03942973
"""  # noqa: E501

# Other invalid cases (these remain unchanged)
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

expected_atoms_force_energy = [
    Atom(atomic_number=16, position=(2.75637007, 2.70758009, 4.82925987)),  # S
    Atom(atomic_number=16, position=(2.84554004, 0.22606000, 3.09241009)),  # S
    Atom(atomic_number=6, position=(1.01258004, 2.49340009, 4.93010997)),  # C
    Atom(atomic_number=6, position=(3.17725992, 1.94638002, 3.18563008)),  # C
    Atom(atomic_number=6, position=(0.00000000, 0.00000000, 0.00000000)),  # C
    Atom(atomic_number=6, position=(2.22169995, -0.06399999, 1.41622007)),  # C
    Atom(atomic_number=6, position=(0.72196996, 0.27478000, 1.31911993)),  # C
    Atom(atomic_number=1, position=(0.73436999, 1.43289995, 4.80946970)),  # H
    Atom(atomic_number=1, position=(0.46778000, 3.08784985, 4.17749977)),  # H
    Atom(atomic_number=1, position=(0.70124000, 2.82674027, 5.93438053)),  # H
    Atom(atomic_number=1, position=(2.59819007, 2.49752998, 2.41988992)),  # H
    Atom(atomic_number=1, position=(4.25116014, 2.06114006, 2.94268989)),  # H
    Atom(atomic_number=1, position=(0.01781000, -1.10166001, -0.05143000)),  # H
    Atom(atomic_number=1, position=(0.50296998, 0.40217999, -0.89547998)),  # H
    Atom(atomic_number=1, position=(-1.05070996, 0.33554000, -0.00612000)),  # H
    Atom(atomic_number=1, position=(2.77840995, 0.53621000, 0.67106002)),  # H
    Atom(atomic_number=1, position=(2.39210987, -1.13260996, 1.18280005)),  # H
    Atom(atomic_number=1, position=(0.59988999, 1.35992002, 1.50082004)),  # H
    Atom(atomic_number=1, position=(0.10381000, -0.26278001, 2.06366992)),  # H
    Atom(atomic_number=8, position=(5.42839002, 0.27989000, 7.42625999)),  # O
    Atom(atomic_number=1, position=(4.56538010, 0.65290999, 7.13586997)),  # H
    Atom(atomic_number=1, position=(5.64198017, 0.76508999, 8.25524997)),  # H
]

expected_gradient = [
    (-1.55959951, 0.22801913, 0.54136018),
    (0.21334177, 1.76227773, -0.33161567),
    (1.31854381, 0.37388875, 0.24927582),
    (0.09910375, -1.52775439, 0.55371287),
    (-0.54047842, 1.02932160, 0.37462899),
    (0.74959058, -0.13849016, 0.66947629),
    (0.88825228, -1.14448169, -0.10462363),
    (-0.11871475, -0.45471041, -0.08255662),
    (-0.24911941, 0.20360289, -0.32966949),
    (-0.20534666, 0.13940533, 0.46775994),
    (-0.29845053, 0.31965051, -0.34195845),
    (0.62592601, -0.53239934, -0.63500463),
    (0.31652062, -0.33261293, 0.62734778),
    (0.14484048, 0.10780390, -0.42454441),
    (-0.55084403, 0.10031390, -0.18450139),
    (0.14551747, 0.35971650, -0.43068238),
    (-0.04825265, -0.56114332, -0.26057253),
    (-0.18014249, 0.45468505, -0.36480505),
    (-0.77321355, -0.35561432, -0.04327458),
    (0.71867266, -0.97644353, -0.67783188),
    (-1.00844588, 0.37843304, -0.30550308),
    (0.31126251, 0.56856677, 1.03720163),
]


expected_molecule = Molecule(
    charge=0,
    multiplicity=1,
    atoms=expected_atoms,
    cell=PeriodicCell(lattice_vectors=expected_cell),
    gradient=None,
)

expected_molecule_with_forces_and_energy = Molecule(
    charge=0,
    multiplicity=1,
    atoms=expected_atoms_force_energy,
    cell=None,
    gradient=expected_gradient,
    energy=-29139.55729148458,
)


def test_molecule_from_extxyz_valid() -> None:
    """
    Test a valid extxyz string (without forces).
    """
    molecule = Molecule.from_extxyz(valid_extxyz)
    assert molecule == expected_molecule, f"Valid case failed:\nGot {molecule}\nExpected {expected_molecule}"


def test_molecule_from_extxyz_valid_with_forces() -> None:
    """
    Test a valid extxyz string that includes forces.
    The forces should be converted to gradients (as the negative of forces).
    """
    molecule = Molecule.from_extxyz(valid_extxyz_with_forces_and_energy)
    assert molecule == expected_molecule_with_forces_and_energy, (
        f"Valid forces case failed:\nGot {molecule}\nExpected {expected_molecule_with_forces_and_energy}"
    )


@mark.parametrize(
    "invalid_extxyz",
    [
        incorrect_num_atoms,
        no_num_atoms,
        not_digit_num_atoms,
        many_num_atoms,
        xyz_style,
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
    Test that invalid extxyz strings raise MoleculeReadError.
    """
    with raises(MoleculeReadError):
        Molecule.from_extxyz(invalid_extxyz)
