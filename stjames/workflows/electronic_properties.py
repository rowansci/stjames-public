from pydantic import NonNegativeInt

from ..base import Base
from ..settings import Settings
from ..types import UUID, ListPerAtom, Matrix3x3, Vector3D
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

    points: list[PropertyCubePoint]


class BondOrder(Base):
    """
    Represents the bond order between two atoms.
    """

    atoms: tuple[NonNegativeInt, NonNegativeInt]
    bond_order: float


class ElectronicPropertiesWorkflow(Workflow):
    """
    Workflow for computing electronic properties!

    Inherited
    :param initial_molecule: molecule of interest

    Config settings:
    :param settings: the level of theory to use
    :param compute_density_cube: whether or not to compute the density on a cube
    :param compute_electrostatic_potential_cube: whether or not to compute the electrostatic potential on a cube
    :param compute_orbital_cubes: which orbitals to compute as cubes. 0 is the lowest energy orbital, 1 is the next lowest in energy, etc.
        For instance, `5` would be water's HOMO, while `6` would be water's LUMO.

    Populated while running:
    :param calculation: the UUID of the calculation
    :param dipole: the dipole moment
    :param quadrupole: the quadrupole moment
    :param mulliken_charges: the Mulliken charges
    :param lowdin_charges: the Lowdin charges
    :param wiberg_bond_orders: the Wiberg bond orders
    :param mayer_bond_orders: the Mayer bond orders
    :param density_cube: the electron density, as a cube
    :param electrostatic_potential_cube: the electrostatic potential, as a cube
    :param orbital_cubes: a dict mapping orbital number to corresponding cube file.
    """

    settings: Settings
    compute_density_cube: bool = True
    compute_electrostatic_potential_cube: bool = True
    compute_orbital_cubes: list[int] = []

    calculation: UUID | None = None

    dipole: Vector3D | None = None
    quadrupole: Matrix3x3 | None = None

    mulliken_charges: ListPerAtom | None = None
    lowdin_charges: ListPerAtom | None = None

    wiberg_bond_orders: list[BondOrder] = []
    mayer_bond_orders: list[BondOrder] = []

    density_cube: PropertyCube | None = None
    electrostatic_potential_cube: PropertyCube | None = None
    orbital_cubes: dict[int, PropertyCube] = {}
