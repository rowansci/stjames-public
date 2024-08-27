from typing import Iterable, Optional, Self, Sequence, TypeAlias

import pydantic
from pydantic import NonNegativeInt, PositiveInt

from .base import Base
from .data import INV_ELEMENT_DICTIONARY, get_atomic_symbol


class MoleculeReadError(RuntimeError):
    pass

PeriodicCell: TypeAlias = tuple[
    tuple[float, float, float],
    tuple[float, float, float],
    tuple[float, float, float],
]


class VibrationalMode(Base):
    frequency: float  # in cm-1
    reduced_mass: float
    force_constant: float
    displacements: list[list[float]]


class Atom(Base):
    atomic_number: NonNegativeInt
    position: list[float]  # in Å

    def __repr__(self) -> str:
        """
        >>> Atom(atomic_number=2, position=[0, 1, 2])
        Atom(2, [0.00000, 1.00000, 2.00000])
        """
        x, y, z = self.position
        return f"Atom({self.atomic_number}, [{x:.5f}, {y:.5f}, {z:.5f}])"

    def __str__(self) -> str:
        """
        >>> str(Atom(atomic_number=2, position=[0, 1, 2]))
        'He    0.0000000000    1.0000000000    2.0000000000'
        """
        x, y, z = self.position
        return f"{self.atomic_symbol:2} {x:15.10f} {y:15.10f} {z:15.10f}"

    @property
    def atomic_symbol(self) -> str:
        """
        >>> Atom(atomic_number=2, position=[0, 1, 2]).atomic_symbol
        'He'
        """
        return get_atomic_symbol(self.atomic_number)

    def edited(self, atomic_number: int | None = None, position: Sequence[float] | None = None) -> Self:
        """
        Create a new Atom with the specified changes.

        >>> a = Atom(atomic_number=2, position=[0, 1, 2])
        >>> a2 = a.edited(3)
        >>> a is a2
        False
        >>> a2
        Atom(3, [0.00000, 1.00000, 2.00000])
        """
        if atomic_number is None:
            atomic_number = self.atomic_number
        if position is None:
            position = list(self.position)

        return self.__class__(atomic_number=atomic_number, position=position)

    @classmethod
    def from_xyz(cls: type[Self], xyz_line: str) -> Self:
        """
        >>> Atom.from_xyz("H 0 0 0")
        Atom(1, [0.00000, 0.00000, 0.00000])
        """
        name, *xyz = xyz_line.split()
        symbol = int(name) if name.isdigit() else INV_ELEMENT_DICTIONARY[name]
        if not len(xyz) == 3:
            raise ValueError("XYZ file should have 3 coordinates per atom")
        return cls(atomic_number=symbol, position=xyz)


class Molecule(Base):
    charge: int
    multiplicity: PositiveInt
    atoms: list[Atom]

    # for periodic boundary conditions
    cell: Optional[PeriodicCell] = None

    energy: Optional[float] = None  # in Hartree
    scf_iterations: Optional[NonNegativeInt] = None
    scf_completed: Optional[bool] = None
    elapsed: Optional[float] = None  # in seconds

    homo_lumo_gap: Optional[float] = None  # in eV

    gradient: Optional[list[list[float]]] = None  # Hartree/Bohr

    mulliken_charges: Optional[list[float]] = None
    mulliken_spin_densities: Optional[list[float]] = None
    dipole: Optional[list[float]] = None  # in Debye

    vibrational_modes: Optional[list[VibrationalMode]] = None

    zero_point_energy: Optional[float] = None
    thermal_energy_corr: Optional[float] = None
    thermal_enthalpy_corr: Optional[float] = None
    thermal_free_energy_corr: Optional[float] = None

    def __len__(self) -> int:
        return len(self.atoms)

    @property
    def coordinates(self) -> list[list[float]]:
        return [a.position for a in self.atoms]

    @property
    def atomic_numbers(self) -> list[NonNegativeInt]:
        return [a.atomic_number for a in self.atoms]

    @property
    def sum_energy_zpe(self) -> Optional[float]:
        if (self.energy is None) or (self.zero_point_energy is None):
            return None
        return self.energy + self.zero_point_energy

    @property
    def sum_energy_thermal_corr(self) -> Optional[float]:
        if (self.energy is None) or (self.thermal_energy_corr is None):
            return None
        return self.energy + self.thermal_energy_corr

    @property
    def sum_energy_enthalpy(self) -> Optional[float]:
        if (self.energy is None) or (self.thermal_enthalpy_corr is None):
            return None
        return self.energy + self.thermal_enthalpy_corr

    @property
    def sum_energy_free_energy(self) -> Optional[float]:
        if (self.energy is None) or (self.thermal_free_energy_corr is None):
            return None
        return self.energy + self.thermal_free_energy_corr

    @pydantic.model_validator(mode="after")
    def check_electron_sanity(self) -> Self:
        num_electrons = sum(self.atomic_numbers) - self.charge
        num_unpaired_electrons = self.multiplicity - 1
        if (num_electrons - num_unpaired_electrons) % 2 != 0:
            raise ValueError(
                f"The combination of {num_electrons} electrons, charge {self.charge}, and multiplicity {self.multiplicity} is impossible. "
                "Double-check the charge and multiplicity values given and verify that they are correct."
            )

        return self

    @pydantic.field_validator("cell")
    @classmethod
    def check_tensor_3D(cls, v: Optional[PeriodicCell]) -> Optional[PeriodicCell]:
        if v is None:
            return v

        if len(v) != 3 or any(len(row) != 3 for row in v):
            raise ValueError("Cell tensor must be a 3x3 list of floats")
        return v

    @classmethod
    def from_xyz(cls: type[Self], xyz: str, charge: int = 0, multiplicity: PositiveInt = 1) -> Self:
        r"""
        Generate a Molecule from an XYZ string.

        Note: only supports single molecule inputs.

        >>> len(Molecule.from_xyz("2\nComment\nH 0 0 0\nH 0 0 1"))
        2
        """
        return cls.from_xyz_lines(xyz.strip().splitlines(), charge=charge, multiplicity=multiplicity)

    @classmethod
    def from_xyz_lines(cls: type[Self], lines: Iterable[str], charge: int = 0, multiplicity: PositiveInt = 1) -> Self:
        lines = list(lines)
        if len(lines[0].split()) == 1:
            natoms = lines[0].strip()
            if not natoms.isdigit() or (int(lines[0]) != len(lines) - 2):
                raise MoleculeReadError(f"First line of XYZ file should be the number of atoms, got: {lines[0]} != {len(lines) - 2}")
            lines = lines[2:]

        try:
            return cls(atoms=[Atom.from_xyz(line) for line in lines], charge=charge, multiplicity=multiplicity)
        except Exception as e:
            raise MoleculeReadError("Error reading molecule from xyz") from e
