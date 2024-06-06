from typing import Optional

from ..base import Base
from .conformer import ConformerSettings
from .workflow import Workflow


class Tautomer(Base):
    energy: float
    weight: Optional[float] = None
    predicted_relative_energy: Optional[float] = None

    # uuids, optionally
    structures: list[Optional[str]] = []


class TautomerWorkflow(Workflow):
    settings: ConformerSettings
    tautomers: list[Tautomer] = []
