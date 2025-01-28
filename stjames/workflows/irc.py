from typing import Self

from pydantic import Field, model_validator

from ..method import XTB_METHODS, Method
from ..mode import Mode
from ..settings import Settings
from ..solvent import Solvent, SolventModel, SolventSettings
from ..types import UUID
from .workflow import Workflow

_sentinel_settings: Settings = object()  # type: ignore [assignment]


class IRCWorkflow(Workflow):
    """
    Workflow for Intrinsic Reaction Coordinate (IRC) calculations.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow

    New:
    :param settings: Settings for running the IRC (only for manual mode)
    :param solvent: Solvent for the calculation (non-Manual mode only)
    :param preopt: whether to optimize the geometry before starting the IRC
    :param final_opt: whether to optimize the final IRC geometry to a minimum

    Results:
    :param starting_TS: optimized TS before the IRC (==initial_molecule if preopt=False)
    :param irc_forward: forward calculations
    :param irc_backward: reverse calculations
    :param opt_forward: optimization steps after the forward IRC
    :param opt_backward: optimization steps after the reverse IRC
    """

    settings: Settings = _sentinel_settings
    solvent: Solvent | None = None

    preopt: bool = False
    final_opt: bool = False

    starting_TS: UUID | None = None

    irc_forward: list[UUID] = Field(default_factory=list)
    irc_backward: list[UUID] = Field(default_factory=list)
    opt_forward: list[UUID] = Field(default_factory=list)
    opt_backward: list[UUID] = Field(default_factory=list)

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        """String representation of the workflow."""
        if self.mode != Mode.MANUAL:
            return f"<{type(self).__name__} {self.mode.name}>"

        return f"<{type(self).__name__} {self.level_of_theory}>"

    @property
    def level_of_theory(self) -> str:
        """Level of theory for the workflow."""
        return self.settings.level_of_theory

    @model_validator(mode="after")
    def validate_mode(self) -> Self:
        """Convert the mode to settings."""
        if self.mode == Mode.MANUAL:
            if self.settings is _sentinel_settings:
                raise ValueError("Settings are required for manual calculations.")
            return self
        elif self.settings is not _sentinel_settings:
            raise ValueError("Cannot specify settings and mode simultaneously.")

        match self.mode:
            case Mode.RAPID:
                method = Method.GFN2_XTB
            case Mode.CAREFUL:
                method = Method.R2SCAN3C
            case Mode.METICULOUS:
                method = Method.WB97X3C
            case _:
                raise ValueError(f"Unsupported mode: {self.mode}")

        model = SolventModel.ALPB if method in XTB_METHODS else SolventModel.CPCM
        solvent_settings = SolventSettings(solvent=self.solvent, model=model) if self.solvent else None

        self.settings = Settings(mode=Mode.RAPID, method=method, solvent_settings=solvent_settings)

        return self
