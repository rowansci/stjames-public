from ..settings import Settings
from ..types import UUID
from .workflow import Workflow


class BasicCalculationWorkflow(Workflow):
    settings: Settings
    engine: str
    calculation_uuid: UUID | None = None
