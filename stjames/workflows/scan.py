from typing import Optional

import numpy as np

from ..base import Base
from ..molecule import Molecule
from ..settings import Settings
from .workflow import Workflow


class ScanPoint(Base):
    index: int
    molecule: Molecule
    energy: Optional[float] = None
    uuid: Optional[str] = None


class ScanSettings(Base):
    type: str  # "bond", "angle", "dihedral" - make Enum later
    atoms: list[int]
    start: float
    stop: float
    num: int

    def vals(self) -> np.ndarray:
        return np.linspace(self.start, self.stop, self.num)

    class Config:
        from_attributes = True


class ScanWorkflow(Workflow):
    scan_settings: ScanSettings
    calc_settings: Settings
    calc_engine: str

    # uuids of scan points
    scan_points: list[str] = []
