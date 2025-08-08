"""Basic calculation workflow."""

from typing import Self

from pydantic import model_validator

from ..engine import Engine
from ..settings import Settings
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
    engine: Engine = None  # type: ignore [assignment]
    calculation_uuid: UUID | None = None

    @model_validator(mode="after")
    def set_engine(self) -> Self:
        """Set the calculation engine."""
        self.engine = self.engine or self.settings.method.default_engine()

        return self
