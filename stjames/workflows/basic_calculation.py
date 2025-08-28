"""Basic calculation workflow."""

from typing import Self

from pydantic import model_validator

from ..base import UniqueList
from ..engine import Engine
from ..settings import Settings
from ..task import Task
from ..types import UUID
from .workflow import MoleculeWorkflow


class BasicCalculationWorkflow(MoleculeWorkflow):
    """
    Workflow for a basic calculation.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow

    New:
    :param settings: Settings for running the calculation
    :param engine: Engine to use
    :param calculation_uuid: UUID of the calculation
    """

    settings: Settings
    tasks: UniqueList[Task] = [Task.ENERGY, Task.CHARGE, Task.DIPOLE]
    calculation_uuid: UUID | None = None

    # DEPRECATED - specify in settings now
    engine: Engine = None  # type: ignore [assignment]

    @model_validator(mode="after")
    def set_engine(self) -> Self:
        """Set the calculation engine."""
        self.engine = self.engine or self.settings.method.default_engine()

        return self
