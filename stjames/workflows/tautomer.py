"""Tautomer prediction workflow."""

from typing import Annotated, Optional

from pydantic import AfterValidator

from ..base import Base, round_float, round_optional_float
from ..mode import Mode
from .workflow import DBCalculation, MoleculeWorkflow


class Tautomer(Base):
    """
    A tautomer.

    :param energy: energy of the tautomer
    :param weight: statistical weight of the tautomer
    :param predicted_relative_energy: relative energy of the tautomer
    :param structures: UUIDs of the structures
    """

    energy: Annotated[float, AfterValidator(round_float(6))]
    weight: Annotated[Optional[float], AfterValidator(round_optional_float(6))] = None
    predicted_relative_energy: Annotated[Optional[float], AfterValidator(round_optional_float(6))] = None

    # UUIDs, optionally
    structures: list[DBCalculation] = []


class TautomerWorkflow(MoleculeWorkflow):
    """
    A workflow to calculate tautomers.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow

    Results:
    :param tautomers: resulting Tautomers
    """

    mode: Mode = Mode.CAREFUL
    tautomers: list[Tautomer] = []
