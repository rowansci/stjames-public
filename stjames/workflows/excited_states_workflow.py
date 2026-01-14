"""Excited State Workflow."""

from typing import Annotated, Literal, Self

from pydantic import BaseModel, Field, PositiveInt, field_validator, model_validator

from stjames.method import DFT_FUNCTIONALS
from stjames.settings import Settings
from stjames.task import Task
from stjames.types import UUID
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
    """

    tda: bool = True
    num_excitations: PositiveInt = 5
    target_root: PositiveInt | None = None

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

        return self


ExcitedStateSettingsUnion = Annotated[TDDFTSettings, Field(discriminator="settings_type")]


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
    calculation_uuid: UUID
