"""Fukui index workflow."""

from typing import Annotated

from pydantic import AfterValidator

from ..base import round_optional_float
from ..method import Method
from ..types import UUID, FloatPerAtom, round_optional_float_per_atom
from .workflow import MoleculeWorkflow


class FukuiIndexWorkflow(MoleculeWorkflow):
    """
    Workflow for calculating Fukui indices.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    Settings:
    :param do_optimization: whether or not to optimize
    :param method: the method
    :param engine_str: the engine

    Results:
    :param optimization: UUID of optimization
    :param global_electrophilicity_index: global electrophilicity index
    :param fukui_positive: Fukui index for positive charges
    :param fukui_negative: Fukui index for negative charges
    :param fukui_zero: Fukui index for zero charges
    """

    do_optimization: bool = True
    optimization: UUID | None = None

    method: Method = Method.GFN1_XTB
    engine_str: str | None = None

    global_electrophilicity_index: Annotated[float | None, AfterValidator(round_optional_float(6))] = None
    fukui_positive: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
    fukui_negative: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
    fukui_zero: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
