from typing import Annotated

from pydantic import AfterValidator

from ..base import round_optional_float
from ..types import UUID
from .workflow import Workflow


class FukuiIndexWorkflow(Workflow):
    # UUID of optimization
    optimization: UUID | None = None

    global_electrophilicity_index: float | None = None
    fukui_positive: Annotated[float | None, AfterValidator(round_optional_float(6))] = None
    fukui_negative: Annotated[float | None, AfterValidator(round_optional_float(6))] = None
    fukui_zero: Annotated[float | None, AfterValidator(round_optional_float(6))] = None
