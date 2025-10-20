"""Fukui index workflow."""

from typing import Annotated, Self

from pydantic import AfterValidator, model_validator

from stjames.method import Method

from ..base import round_optional_float
from ..engine import Engine
from ..settings import Settings
from ..types import UUID, FloatPerAtom, round_optional_float_per_atom
from .workflow import MoleculeWorkflow


class FukuiIndexWorkflow(MoleculeWorkflow):
    """
    Workflow for calculating Fukui indices.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    Settings:
    :param opt_settings: if given, the settings for optimization. if none, no optimization will be conducted.
    :param opt_engine: the engine for optimization [uses opt_settings.method.default_engine if not set]
    :param fukui_settings: the settings for Fukui index calculations.
    :param fukui_engine: the engine for Fukui index calculations [uses fukui_settings.method.default_engine if not set]

    Results:
    :param optimization: UUID of optimization
    :param global_electrophilicity_index: global electrophilicity index
    :param fukui_positive: Fukui index for positive charges
    :param fukui_negative: Fukui index for negative charges
    :param fukui_zero: Fukui index for zero charges
    """

    opt_settings: Settings | None = None
    opt_engine: Engine | None = None

    fukui_settings: Settings = Settings(method=Method.GFN1_XTB)
    fukui_engine: Engine = None  # type: ignore [assignment]

    optimization: UUID | None = None

    global_electrophilicity_index: Annotated[float | None, AfterValidator(round_optional_float(6))] = None
    fukui_positive: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
    fukui_negative: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None
    fukui_zero: Annotated[FloatPerAtom | None, AfterValidator(round_optional_float_per_atom(6))] = None

    @model_validator(mode="after")
    def set_engines(self) -> Self:
        """Set the engines for optimization and Fukui index calculations."""
        if self.opt_settings is not None and self.opt_engine is None:
            self.opt_engine = self.opt_settings.method.default_engine()
        self.fukui_engine = self.fukui_engine or self.fukui_settings.method.default_engine()

        return self
