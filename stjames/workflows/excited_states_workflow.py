"""Excited State Workflow."""

from typing import Annotated, Literal, Self

from pydantic import AfterValidator, BaseModel, Field, PositiveFloat, PositiveInt, field_validator, model_validator

from stjames.base import LowercaseStrEnum, round_optional_float
from stjames.method import DFT_FUNCTIONALS, RANGE_SEPARATED_FUNCTIONALS
from stjames.molecule import Molecule
from stjames.settings import Settings
from stjames.task import Task
from stjames.types import round_list
from stjames.workflows.workflow import MoleculeWorkflow


class ExcitedStateSettings(BaseModel):
    """
    Settings for excited state calculations.

    :param settings: calculation settings (used for method, basis set, etc.)
    :param tasks: Tasks to perform
    """

    settings: Settings
    tasks: set[Literal[Task.ENERGY, Task.GRADIENT, Task.OPTIMIZE, Task.FREQUENCIES]] = {Task.ENERGY}

    @field_validator("tasks")
    @classmethod
    def validate_tasks(cls, v: set[Task]) -> set[Task]:
        """Remove GRADIENT if OPTIMIZE is also present."""
        if Task.GRADIENT in v and Task.OPTIMIZE in v:
            v.remove(Task.GRADIENT)

        return v


class OmegaTuning(LowercaseStrEnum):
    """Options for omega (range-separation parameter) tuning"""

    KOOPMANS = "koopmans"  # Baer et al. doi.org/10.1146/annurev.physchem.012809.103321


class TDDFTSettings(ExcitedStateSettings):
    """
    Settings for TDDFT calculations.

    Inherited:
    :param settings: calculation settings (used for method, basis set, etc.)
    :param tasks: Tasks to perform

    New:
    :param tda: use Tamm-Dancoff approximation
    :param num_excitations: number of excitations to calculate
    :param target_root: root to target (for gradient/optimization)
    :param omega: range-separation parameter (Bohr⁻¹) or method to tune it (optional)
    """

    tda: bool = True
    num_excitations: PositiveInt = 5
    target_root: PositiveInt | None = None
    omega: PositiveFloat | OmegaTuning | None = None

    settings_type: Literal["TDDFTSettings"] = "TDDFTSettings"

    @field_validator("settings")
    def validate_settings(cls, settings: Settings) -> Settings:
        """Validate setup"""
        if settings.method not in DFT_FUNCTIONALS:
            raise ValueError("TDDFT can only be run with DFT functionals.")

        return settings

    @model_validator(mode="after")
    def validate_setup(self) -> Self:
        if self.tasks & {Task.GRADIENT, Task.OPTIMIZE, Task.FREQUENCIES}:
            if self.target_root is None:
                raise ValueError("target_root must be specified when GRADIENT or OPTIMIZE task is selected.")
            if self.target_root > self.num_excitations:
                raise ValueError("target_root cannot be greater than num_excitations.")
        elif self.target_root is not None:
            raise ValueError("target_root requires GRADIENT, OPTIMIZE, or FREQUENCIES task to be selected.")

        if self.omega is not None and self.settings.method not in RANGE_SEPARATED_FUNCTIONALS:
            functionals = "\n    ".join(RANGE_SEPARATED_FUNCTIONALS)
            raise ValueError(f"Omega may only be specified for range-separated DFT functionals:\n    {functionals}.")

        return self


class ExcitedStateResult(BaseModel):
    """
    Results of an excited state calculation.

    :param molecule: Molecule containing ground-state energy and excited state gradient/frequencies where applicable
    """

    molecule: Molecule


class TDDFTResult(ExcitedStateResult):
    """
    Results of a TDDFT calculation.

    Inherited:
    :param molecule: Molecule containing ground-state energy and excited state gradient/frequencies where applicable

    New:
    :param omega: range-separation parameter (Bohr⁻¹), if omega tuning was performed or different than the base functional
    :param excitation_energies: excitation energies in Hartree
    :param oscillator_strengths: oscillator strengths (length where applicable)
    """

    omega: Annotated[float | None, AfterValidator(round_optional_float(5))] = None
    excitation_energies: Annotated[list[float], AfterValidator(round_list(3))]
    oscillator_strengths: Annotated[list[float], AfterValidator(round_list(3))]

    result_type: Literal["TDDFTResult"] = "TDDFTResult"

    @model_validator(mode="after")
    def validate_results(self) -> Self:
        if len(self.excitation_energies) != len(self.oscillator_strengths):
            raise ValueError("excitation_energies and oscillator_strengths must have the same length.")

        return self


ExcitedStateSettingsUnion = Annotated[TDDFTSettings, Field(discriminator="settings_type")]
ExcitedStateResultUnion = Annotated[TDDFTResult, Field(discriminator="result_type")]


class ExcitedStatesWorkflow(MoleculeWorkflow):
    """
    Workflow for an excited state calculation.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (unused)

    New:
    :param excited_state_settings: settings for running the excited state calculation
    :param results: results of the excited state calculation
    """

    excited_state_settings: ExcitedStateSettingsUnion
    results: list[ExcitedStateResultUnion] | None = None
