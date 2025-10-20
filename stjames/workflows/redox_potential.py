"""Redox potential workflow."""

from typing import Annotated, Any, TypeVar

from pydantic import AfterValidator, ValidationInfo, field_validator, model_validator

from ..base import round_optional_float
from ..mode import Mode
from ..solvent import Solvent
from ..types import UUID
from .multistage_opt import MultiStageOptMixin
from .workflow import MoleculeWorkflow

_T = TypeVar("_T")


class RedoxPotentialWorkflow(MoleculeWorkflow, MultiStageOptMixin):
    """
    Workflow for computing spin states of molecules.

    Uses the modes from MultiStageOptSettings.

    Inherited
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow
    :param multistage_opt_settings: set by mode unless mode=MANUAL (ignores additional settings if set)
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param constraints: constraints to add
    :param transition_state: whether this is a transition state
    :param frequencies: whether to calculate frequencies

    Overridden:
    :param mso_mode: Mode for MultiStageOptSettings
    :param solvent: solvent to use for optimization

    New:
    :param reduction: whether or not to calculate the reduction half-reaction
    :param oxidation: whether or not to calculate the oxidation half-reaction
    :param neutral_molecule: UUID of the calculation for the neutral molecule
    :param anion_molecule: UUID of the calculation for the anion molecule
    :param cation_molecule: UUID of the calculation for the cation molecule
    :param reduction_potential: the final potential, in V
    :param oxidation_potential: the final potential, in V

    Legacy:
    :param redox_type: one of "reduction" or "oxidation"
    :param redox_potential: the corresponding potential, in V
    """

    solvent: Solvent = Solvent.ACETONITRILE

    reduction: bool = True
    oxidation: bool = True

    # legacy values - remove in future release!
    redox_type: str | None = None
    redox_potential: float | None = None

    # UUIDs
    neutral_molecule: UUID | None = None
    anion_molecule: UUID | None = None
    cation_molecule: UUID | None = None

    reduction_potential: Annotated[float | None, AfterValidator(round_optional_float(6))] = None
    oxidation_potential: Annotated[float | None, AfterValidator(round_optional_float(6))] = None

    @field_validator("solvent", mode="before")
    @classmethod
    def only_mecn_please(cls, val: Solvent | None) -> Solvent:
        """Only MeCN please!"""
        if val != Solvent.ACETONITRILE:
            raise ValueError("Only acetonitrile permitted!")

        return val

    @field_validator("constraints", "transition_state")
    @classmethod
    def turned_off(cls, value: _T, info: ValidationInfo) -> _T:
        if value:
            raise ValueError(f"{info.field_name} not supported in redox potential workflows.")

        return value

    @model_validator(mode="before")
    @classmethod
    def set_mode_and_mso_mode(cls, values: dict[str, Any]) -> dict[str, Any]:
        """Set the MultiStageOptSettings mode to match current redox potential mode, and select mode if `Auto`."""
        if ("mode" not in values) or (values["mode"] == Mode.AUTO):
            values["mode"] = Mode.RAPID

        values["mso_mode"] = values["mode"]
        return values

    def model_post_init(self, __context: Any) -> None:
        """Keep back-compatible with old schema."""
        if self.redox_type == "oxidation":
            self.oxidation = True
            self.reduction = False
            if self.oxidation_potential is None:
                self.oxidation_potential = self.redox_potential
        elif self.redox_type == "reduction":
            self.oxidation = False
            self.reduction = True
            if self.reduction_potential is None:
                self.reduction_potential = self.redox_potential
