from typing import Annotated, Optional

from pydantic import AfterValidator

from ..base import Base, round_float, round_optional_float
from ..mode import Mode
from .workflow import DBCalculation, Workflow


class Tautomer(Base):
    energy: Annotated[float, AfterValidator(round_float(6))]
    weight: Annotated[Optional[float], AfterValidator(round_optional_float(6))] = None
    predicted_relative_energy: Annotated[Optional[float], AfterValidator(round_optional_float(6))] = None

    # UUIDs, optionally
    structures: list[DBCalculation] = []


class TautomerWorkflow(Workflow):
    mode: Mode = Mode.CAREFUL
    tautomers: list[Tautomer] = []
