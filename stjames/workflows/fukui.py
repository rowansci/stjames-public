from typing import Optional

from .workflow import Workflow


class FukuiIndexWorkflow(Workflow):
    # uuid of optimization
    optimization: Optional[str] = None

    global_electrophilicity_index: float
    fukui_positive: list[float]
    fukui_negative: list[float]
    fukui_zero: list[float]
