"""Scan workflow."""

from typing import Annotated

import numpy as np
from numpy.typing import NDArray
from pydantic import AfterValidator, field_validator

from ..base import Base, round_optional_float
from ..molecule import Molecule
from ..settings import Settings
from ..types import UUID
from .workflow import MoleculeWorkflow


class ScanPoint(Base):
    """
    A point in a scan.

    :param index: index of the point
    :param index_2d: index of the point in 2nd dimension, if applicable
    :param molecule: Molecule at the point
    :param energy: energy of the point
    :param uuid: UUID of the calculation
    """

    index: int
    index_2d: int = 0
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
    :param scan_settings: what coordinate(s) to scan; if more than one, all will be performed simultaneously and should have the same number of steps
    :param scan_settings_2d: what additional coordinate(s) to scan; makes a grid with `scan_settings`
    :param wavefront propagation: whether to use wavefront propagation (10.1063/5.0009232) for more expensive but smoother scans
    :param calc_settings: settings for the calculation
    :param calc_engine: engine to use for the calculation
    :param scan_points: points along the scan
    """

    scan_settings: ScanSettings | list[ScanSettings]
    scan_settings_2d: ScanSettings | list[ScanSettings] = []
    calc_settings: Settings
    calc_engine: str

    wavefront_propagation: bool = True

    # UUIDs of scan points
    scan_points: list[UUID | None] = []

    @field_validator("scan_settings", "scan_settings_2d", mode="after")
    @classmethod
    def validate_scan_settings(cls, val: ScanSettings | list[ScanSettings]) -> list[ScanSettings]:
        """Ensure that scan_settings is a list, and that every list item has the same number of steps."""
        if isinstance(val, ScanSettings):
            val = [val]

        for ss in val:
            if ss.num != val[0].num:
                raise ValueError("Concerted scan settings must have same number of steps!")

        return val
