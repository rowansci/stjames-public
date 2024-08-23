from typing import Optional

from .workflow import Workflow


class ADMETWorkflow(Workflow):
    properties: Optional[dict[str, float | int]] = None
