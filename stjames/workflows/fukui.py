"""Fukui index workflow."""

from typing import Annotated

from pydantic import AfterValidator

from ..base import round_optional_float
from ..types import UUID, FloatPerAtom, round_optional_float_per_atom
from .workflow import MoleculeWorkflow


class FukuiIndexWorkflow(MoleculeWorkflow):
    """
    Workflow for calculating Fukui indices.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    Results:
    :param optimization: UUID of optimization
    :param global_electrophilicity_index: global electrophilicity index
    :param fukui_positive: Fukui index for positive charges
    :param fukui_negative: Fukui index for negative charges
    :param fukui_zero: Fukui index for zero charges
    """

    # UUID of optimization
    optimization: UUID | None = None

    global_electrophilicity_index: Annotated[float | None, AfterValidator(round_optional_float(6))] = None
    fukui_positive: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
    fukui_negative: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
    fukui_zero: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
