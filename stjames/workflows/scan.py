"""Scan workflow."""

from typing import Annotated

import numpy as np
from numpy.typing import NDArray
from pydantic import AfterValidator

from ..base import Base, round_optional_float
from ..molecule import Molecule
from ..settings import Settings
from ..types import UUID
from .workflow import MoleculeWorkflow


class ScanPoint(Base):
    """
    A point in a scan.

    :param index: index of the point
    :param molecule: Molecule at the point
    :param energy: energy of the point
    :param uuid: UUID of the calculation
    """

    index: int
    molecule: Molecule
    energy: Annotated[float | None, AfterValidator(round_optional_float(6))] = None
    uuid: UUID | None = None


class ScanSettings(Base):
    """
    Settings for a scan.

    :param type: type of scan (bond, angle, dihedral)
    :param atoms: indices of atoms to scan
    :param start: start of scan
    :param stop: end of scan
    :param num: number of points in the scan
    """

    type: str  # "bond", "angle", "dihedral" - make Enum later
    atoms: list[int]
    start: float
    stop: float
    num: int

    def vals(self) -> NDArray[np.float64]:
        return np.linspace(self.start, self.stop, self.num)  # type: ignore [return-value, unused-ignore]

    class Config:
        from_attributes = True


class ScanWorkflow(MoleculeWorkflow):
    """
    Workflow for scanning a bond, angle, or dihedral.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param scan_settings: information about what coordinate to scan
    :param calc_settings: settings for the calculation
    :param calc_engine: engine to use for the calculation
    :param scan_points: points along the scan
    """

    scan_settings: ScanSettings
    calc_settings: Settings
    calc_engine: str

    # UUIDs of scan points
    scan_points: list[UUID | None] = []
