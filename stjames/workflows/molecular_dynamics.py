from typing import Self

from pydantic import PositiveFloat, PositiveInt, computed_field, model_validator

from ..base import Base, LowercaseStrEnum
from ..constraint import PairwiseHarmonicConstraint, SphericalHarmonicConstraint
from ..settings import Settings
from ..types import UUID
from .workflow import Workflow


class MolecularDynamicsInitialization(LowercaseStrEnum):
    RANDOM = "random"
    QUASICLASSICAL = "quasiclassical"
    READ = "read"


class ThermodynamicEnsemble(LowercaseStrEnum):
    NPT = "npt"
    NVT = "nvt"
    NVE = "nve"


class Frame(Base):
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


class MolecularDynamicsWorkflow(Workflow):
    settings: MolecularDynamicsSettings
    calc_settings: Settings
    calc_engine: str | None = None

    # UUIDs of scan points
    frames: list[Frame] = []
