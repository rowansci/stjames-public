from typing import Annotated

from pydantic import AfterValidator

from ..base import round_optional_float
from ..types import UUID, FloatPerAtom, round_optional_float_per_atom
from .workflow import Workflow


class FukuiIndexWorkflow(Workflow):
    # UUID of optimization
    optimization: UUID | None = None

    global_electrophilicity_index: Annotated[float | None, AfterValidator(round_optional_float(6))] = None
    fukui_positive: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
    fukui_negative: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
    fukui_zero: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
