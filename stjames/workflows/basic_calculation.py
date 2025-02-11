"""Basic calculation workflow."""

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
    engine: str
    calculation_uuid: UUID | None = None
