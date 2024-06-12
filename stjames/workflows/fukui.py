from typing import Optional

from .workflow import Workflow


class FukuiIndexWorkflow(Workflow):
    # uuid of optimization
    optimization: Optional[str] = None

    global_electrophilicity_index: Optional[float] = None
    fukui_positive: Optional[list[float]] = None
    fukui_negative: Optional[list[float]] = None
    fukui_zero: Optional[list[float]] = None
