from pytest import fixture

from stjames import Atom, ElectronicPropertiesWorkflow, Molecule, Settings


@fixture
def Mn() -> Molecule:
    return Molecule(charge=0, multiplicity=2, atoms=[Atom(atomic_number=25, position=[0, 0, 0])])


@fixture
def Fe() -> Molecule:
    return Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=26, position=[0, 0, 0])])


def test_electronic_properties(Mn: Molecule, Fe: Molecule) -> None:
    epw1 = ElectronicPropertiesWorkflow(
        initial_molecule=Mn,
        settings=Settings(
            functional="PBE",
            basis="def2-SVP",
        ),
    )

    assert epw1.compute_density_cube is True
    assert epw1.compute_electrostatic_potential_cube is True
    assert epw1.compute_num_occupied_orbitals == 1
    assert epw1.compute_num_virtual_orbitals == 1
    assert epw1.calc_uuid is None
    assert epw1.dipole is None
    assert epw1.quadrupole is None
    assert epw1.lowdin_charges is None
    assert epw1.mulliken_charges is None
    assert epw1.wiberg_bond_orders == []
    assert epw1.mayer_bond_orders == []
    assert epw1.density_cube is None
    assert epw1.density_cube_alpha is None
    assert epw1.density_cube_beta is None
    assert epw1.density_cube_difference is None
    assert epw1.electrostatic_potential_cube is None
    assert epw1.molecular_orbitals == {}
    assert epw1.molecular_orbitals_alpha == {}
    assert epw1.molecular_orbitals_beta == {}

    epw2 = ElectronicPropertiesWorkflow(
        initial_molecule=Fe,
        settings=Settings(
            functional="PBEh",
            basis="def2-SVP",
        ),
        compute_density_cube=False,
        compute_electrostatic_potential_cube=False,
        compute_num_occupied_orbitals=2,
        compute_num_virtual_orbitals=2,
    )

    assert epw2.compute_density_cube is False
    assert epw2.compute_electrostatic_potential_cube is False
    assert epw2.compute_num_occupied_orbitals == 2
    assert epw2.compute_num_virtual_orbitals == 2
    assert epw2.calc_uuid is None
    assert epw2.dipole is None
    assert epw2.quadrupole is None
    assert epw2.lowdin_charges is None
    assert epw2.mulliken_charges is None
    assert epw2.wiberg_bond_orders == []
    assert epw2.mayer_bond_orders == []
    assert epw2.density_cube is None
    assert epw2.density_cube_alpha is None
    assert epw2.density_cube_beta is None
    assert epw2.density_cube_difference is None
    assert epw2.electrostatic_potential_cube is None
    assert epw2.molecular_orbitals == {}
    assert epw2.molecular_orbitals_alpha == {}
    assert epw2.molecular_orbitals_beta == {}
