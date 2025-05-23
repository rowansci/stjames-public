"""Molecular dynamics workflow."""

from typing import Self

from pydantic import PositiveFloat, PositiveInt, computed_field, model_validator

from ..base import Base, LowercaseStrEnum
from ..constraint import PairwiseHarmonicConstraint, SphericalHarmonicConstraint
from ..settings import Settings
from ..types import UUID
from .workflow import MoleculeWorkflow


class MolecularDynamicsInitialization(LowercaseStrEnum):
    """Initialization method for molecular dynamics."""

    RANDOM = "random"
    QUASICLASSICAL = "quasiclassical"
    READ = "read"


class ThermodynamicEnsemble(LowercaseStrEnum):
    """Thermodynamic ensemble for molecular dynamics."""

    NPT = "npt"
    NVT = "nvt"
    NVE = "nve"


class Frame(Base):
    """A frame in a molecular dynamics simulation."""

    index: int  # what number frame this is within the MD simulation

    calculation_uuid: UUID | None = None  # UUID of calculation

    pressure: float
    temperature: float
    volume: float
    potential_energy: float  # kcal/mol
    kinetic_energy: float  # kcal/mol

    @computed_field  # type: ignore[misc, prop-decorator, unused-ignore]
    @property
    def energy(self) -> float:
        return self.potential_energy + self.kinetic_energy


class MolecularDynamicsSettings(Base):
    """Settings for a molecular dynamics simulation."""

    ensemble: ThermodynamicEnsemble = ThermodynamicEnsemble.NVT
    initialization: MolecularDynamicsInitialization = MolecularDynamicsInitialization.RANDOM

    timestep: PositiveFloat = 1.0  # fs
    num_steps: PositiveInt = 500
    save_interval: PositiveInt = 10

    confining_constraint: SphericalHarmonicConstraint | None = None

    temperature: PositiveFloat = 300  # K
    pressure: PositiveFloat | None = None  # atm

    langevin_thermostat_timescale: PositiveFloat = 100  # fs
    berendsen_barostat_timescale: PositiveFloat = 1000  # fs

    constraints: list[PairwiseHarmonicConstraint] = []

    @model_validator(mode="after")
    def validate_ensemble_settings(self) -> Self:
        """Check that NVT ensemble always has temperature defined, and that NPT has temp and pressure defined."""
        if self.ensemble == ThermodynamicEnsemble.NVT and self.temperature is None:
            raise ValueError("NVT ensemble must have a temperature defined")
        elif self.ensemble == ThermodynamicEnsemble.NPT and (self.temperature is None or self.pressure is None):
            raise ValueError("NPT ensemble must have both temperature and pressure defined")
        return self


class MolecularDynamicsWorkflow(MoleculeWorkflow):
    """
    Workflow for running a molecular dynamics simulation.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param settings: settings for the molecular dynamics simulation
    :param calc_settings: settings for the gradient calculation
    :param calc_engine: engine to use for the gradient calculation

    Results:
    :param frames: Frames from the MD
    """

    settings: MolecularDynamicsSettings
    calc_settings: Settings
    calc_engine: str | None = None

    frames: list[Frame] = []
