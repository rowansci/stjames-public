from stjames import Atom, Molecule, VibrationalMode


def test_vibrational_mode_rounding() -> None:
    test_vectors = [(1.23456789345346, 2.34567890346346, 3.45678901346346436), (4.567890126343636, 5.678901236346436, 6.7890123436436436)]
    mode = VibrationalMode(frequency=1.1234567, reduced_mass=1.1234567, force_constant=1.1234567, displacements=test_vectors)

    expected_vectors = [(1.234568, 2.345679, 3.456789), (4.567890, 5.678901, 6.789012)]
    assert mode.displacements == expected_vectors
    assert mode.frequency == 1.123
    assert mode.reduced_mass == 1.123
    assert mode.force_constant == 1.123


def test_atom_rounding() -> None:
    atom = Atom(atomic_number=2, position=[0.1111111112222222, 1.1111111112222222, 2.1111111112222222])

    rounded_position = (0.11111111, 1.11111111, 2.11111111)

    assert atom.position == rounded_position


def test_molecule_rounding() -> None:
    mol = Molecule(
        charge=0,
        multiplicity=1,
        atoms=[Atom(atomic_number=2, position=[0.1111111112222222, 1.1111111112222222, 2.1111111112222222])],
        energy=1.23456789345346,
        zero_point_energy=2.34567890346346,
        thermal_energy_corr=3.45678901346346436,
        elapsed=4.567890126343636,
        thermal_enthalpy_corr=3.45678901346346436,
        thermal_free_energy_corr=3.45678901346346436,
        homo_lumo_gap=5.12345678,
        dipole=(1.23456789345346, 2.34567890346346, 3.45678901346346436),
        stress=(
            (1.23456789345346, 2.34567890346346, 3.45678901346346436),
            (4.567890126343636, 5.678901236346436, 6.7890123436436436),
            (4.567890126343636, 5.678901236346436, 6.7890123436436436),
        ),
    )

    assert mol.atoms[0].position == (0.11111111, 1.11111111, 2.11111111)
    assert mol.energy == 1.234568
    assert mol.zero_point_energy == 2.345679
    assert mol.thermal_energy_corr == 3.456789
    assert mol.elapsed == 4.568
    assert mol.thermal_enthalpy_corr == 3.456789
    assert mol.thermal_free_energy_corr == 3.456789
    assert mol.homo_lumo_gap == 5.123457
    assert mol.dipole == (1.234568, 2.345679, 3.456789)
    assert mol.stress == ((1.234568, 2.345679, 3.456789), (4.567890, 5.678901, 6.789012), (4.567890, 5.678901, 6.789012))
