from typing import Annotated

from pydantic import AfterValidator, NonNegativeFloat, NonNegativeInt

from ..base import Base, round_float
from ..settings import Settings
from ..types import UUID, FloatPerAtom, Matrix3x3, Vector3D, round_optional_float_per_atom, round_optional_matrix3x3, round_optional_vector3d
from .workflow import Workflow


class PropertyCubePoint(Base):
    """A point in a cube file, all values rounded to 6 decimal places."""

    x: Annotated[float, AfterValidator(round_float(3))]
    y: Annotated[float, AfterValidator(round_float(3))]
    z: Annotated[float, AfterValidator(round_float(3))]
    val: Annotated[float, AfterValidator(round_float(6))]


class PropertyCube(Base):
    """Represents a "cubefile" of some property."""

    cube_points: list[PropertyCubePoint]


class MolecularOrbitalCube(PropertyCube):
    """
    Cube of a molecular orbital.

    Inherits `cube_data`.
    """

    occupation: NonNegativeInt
    energy: Annotated[float, AfterValidator(round_float(6))]


class ElectronicPropertiesWorkflow(Workflow):
    """
    Workflow for computing electronic properties.

    Inherited
    :param initial_molecule: Molecule of interest

    Config settings:
    :param settings: settings for the calculation
    :param compute_density_cube: whether to compute the density cube
    :param compute_electrostatic_potential_cube: whether to compute the electrostatic potential cube
    :param compute_num_occupied_orbitals: number of occupied orbitals to save
    :param compute_num_virtual_orbitals: number of virtual orbitals to save

    Populated while running:
    :param calc_uuid: UUID of the calculation
    :param dipole: dipole moment
    :param quadrupole: quadrupole moment
    :param lowdin_charges: Löwdin charges
    :param mulliken_charges: Mulliken charges
    :param wiberg_bond_orders: Wiberg bond orders (`atom1`, `atom2`, `order`)
    :param mayer_bond_orders: Mayer bond orders (`atom1`, `atom2`, `order`)

    :param density_cube: electron density, as a cube
    :param density_cube_alpha: α electron density, as a cube
    :param density_cube_beta: β electron density, as a cube
    :param density_cube_difference: difference spin densities, as a cube

    :param electrostatic_potential_cube: electrostatic potential, as a cube

    :param molecular_orbitals: MOs, key is absolute orbital index (for closed-shell species (RHF))
    :param molecular_orbitals_alpha: α MOs, key is absolute orbital index (for open-shell species (UHF/ROHF))
    :param molecular_orbitals_beta: β MOs, key is absolute orbital index (for open-shell species (UHF/ROHF))
    """

    # Config settings
    settings: Settings
    compute_density_cube: bool = True
    compute_electrostatic_potential_cube: bool = True
    compute_num_occupied_orbitals: NonNegativeInt = 1
    compute_num_virtual_orbitals: NonNegativeInt = 1

    # Results
    calc_uuid: UUID | None = None

    dipole: Annotated[Vector3D | None, AfterValidator(round_optional_vector3d(6))] = None
    quadrupole: Annotated[Matrix3x3 | None, AfterValidator(round_optional_matrix3x3(6))] = None

    mulliken_charges: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
    lowdin_charges: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None

    wiberg_bond_orders: list[tuple[NonNegativeInt, NonNegativeInt, NonNegativeFloat]] = []
    mayer_bond_orders: list[tuple[NonNegativeInt, NonNegativeInt, NonNegativeFloat]] = []

    density_cube: PropertyCube | None = None
    density_cube_alpha: PropertyCube | None = None
    density_cube_beta: PropertyCube | None = None
    density_cube_difference: PropertyCube | None = None

    electrostatic_potential_cube: PropertyCube | None = None

    molecular_orbitals: dict[NonNegativeInt, MolecularOrbitalCube] = {}
    molecular_orbitals_alpha: dict[NonNegativeInt, MolecularOrbitalCube] = {}
    molecular_orbitals_beta: dict[NonNegativeInt, MolecularOrbitalCube] = {}
