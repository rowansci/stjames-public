from ..types import UUID
from .workflow import Workflow


class FukuiIndexWorkflow(Workflow):
    # uuid of optimization
    optimization: UUID | None = None

    global_electrophilicity_index: float | None = None
    fukui_positive: list[float] | None = None
    fukui_negative: list[float] | None = None
    fukui_zero: list[float] | None = None
