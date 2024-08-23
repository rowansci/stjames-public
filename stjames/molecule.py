from typing import Optional, Self

import pydantic
from pydantic import NonNegativeInt, PositiveInt

from .base import Base


class VibrationalMode(Base):
    frequency: float  # in cm-1
    reduced_mass: float
    force_constant: float
    displacements: list[list[float]]


class Atom(Base):
    atomic_number: NonNegativeInt
    position: list[float]  # in Å


class Molecule(Base):
    charge: int
    multiplicity: PositiveInt
    atoms: list[Atom]

    # for periodic boundary conditions
    cell: Optional[list[list[float]]] = None

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
    def check_tensor(cls, v: list[list[float]]) -> list[list[float]]:
        if len(v) != 3 or any(len(row) != 3 for row in v):
            raise ValueError("Cell tensor must be a 3x3 list of floats")
        return v
