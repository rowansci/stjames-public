from typing import Optional

from ..settings import Settings
from .workflow import Workflow


class BasicCalculationWorkflow(Workflow):
    settings: Settings
    engine: str
    calculation_uuid: Optional[str] = None
