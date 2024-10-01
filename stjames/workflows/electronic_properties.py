from pydantic import NonNegativeFloat, NonNegativeInt

from ..base import Base
from ..settings import Settings
from ..types import UUID, FloatPerAtom, Matrix3x3, Vector3D
from .workflow import Workflow


class PropertyCubePoint(Base):
    x: float
    y: float
    z: float
    val: float


class PropertyCube(Base):
    """
    Represents a "cubefile" of some property.
    """

    cube_data: list[PropertyCubePoint]


class MolecularOrbitalCube(PropertyCube):
    """
    Inherits `cube_data`.
    """

    occupation: NonNegativeInt
    energy: float


class ElectronicPropertiesWorkflow(Workflow):
    """
    Workflow for computing electronic properties!

    Inherited
    :param initial_molecule: molecule of interest

    Config settings:
    :param settings: the level of theory to use
    :param compute_density_cube: whether or not to compute the density on a cube
    :param compute_electrostatic_potential_cube: whether or not to compute the electrostatic potential on a cube
    :param compute_num_occupied_orbitals: how many occupied orbitals to save
    :param compute_num_virtual_orbitals: how many virtual orbitals to save

    Populated while running:
    :param calculation: the UUID of the calculation
    :param dipole: the dipole moment
    :param quadrupole: the quadrupole moment
    :param mulliken_charges: the Mulliken charges
    :param lowdin_charges: the Lowdin charges
    :param wiberg_bond_orders: the Wiberg bond orders (`atom1`, `atom2`, `order`)
    :param mayer_bond_orders: the Mayer bond orders (`atom1`, `atom2`, `order`)
    :param density_cube: the electron density, as a cube
    :param electrostatic_potential_cube: the electrostatic potential, as a cube
    :param molecular_orbitals_alpha: for open-shell species (UHF/ROHF), a dict containing the alpha MOs
        (The key is the absolute orbital index.)
    :param molecular_orbitals_beta: for open-shell species (UHF/ROHF), a dict containing the beta MOs
        (The key is the absolute orbital index.)
    :param molecular_orbitals: for closed-shell species (RHF), a dict containing the MOs
        (The key is the absolute orbital index.)
    """

    settings: Settings
    compute_density_cube: bool = True
    compute_electrostatic_potential_cube: bool = True
    compute_num_occupied_orbitals: NonNegativeInt = 1
    compute_num_virtual_orbitals: NonNegativeInt = 1

    calculation: UUID | None = None

    dipole: Vector3D | None = None
    quadrupole: Matrix3x3 | None = None

    mulliken_charges: FloatPerAtom | None = None
    lowdin_charges: FloatPerAtom | None = None

    wiberg_bond_orders: list[tuple[NonNegativeInt, NonNegativeInt, NonNegativeFloat]] = []
    mayer_bond_orders: list[tuple[NonNegativeInt, NonNegativeInt, NonNegativeFloat]] = []

    density_cube: PropertyCube | None = None
    electrostatic_potential_cube: PropertyCube | None = None
    molecular_orbitals_alpha: dict[NonNegativeInt, MolecularOrbitalCube] = {}
    molecular_orbitals_beta: dict[NonNegativeInt, MolecularOrbitalCube] = {}
    molecular_orbitals: dict[NonNegativeInt, MolecularOrbitalCube] = {}
