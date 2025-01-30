from typing import Annotated

import numpy as np
from numpy.typing import NDArray
from pydantic import AfterValidator

from ..base import Base, round_optional_float
from ..molecule import Molecule
from ..settings import Settings
from ..types import UUID
from .workflow import Workflow


class ScanPoint(Base):
    index: int
    molecule: Molecule
    energy: Annotated[float | None, AfterValidator(round_optional_float(6))] = None
    uuid: UUID | None = None


class ScanSettings(Base):
    type: str  # "bond", "angle", "dihedral" - make Enum later
    atoms: list[int]
    start: float
    stop: float
    num: int

    def vals(self) -> NDArray[np.float64]:
        return np.linspace(self.start, self.stop, self.num)  # type: ignore [return-value, unused-ignore]

    class Config:
        from_attributes = True


class ScanWorkflow(Workflow):
    scan_settings: ScanSettings
    calc_settings: Settings
    calc_engine: str

    # UUIDs of scan points
    scan_points: list[UUID | None] = []
