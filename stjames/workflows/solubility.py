"""Solubility prediction workflow."""

from typing import Annotated, Self

from pydantic import AfterValidator, BaseModel, model_validator

from ..base import LowercaseStrEnum
from ..types import round_list
from .workflow import SMILESWorkflow


class SolubilityMethod(LowercaseStrEnum):
    FASTSOLV = "fastsolv"
    KINGFISHER = "kingfisher"
    ESOL = "esol"


class SolubilityResult(BaseModel):
    """
    Solubility result.

    :param solubilities: solubilities in log(S)/L
    :param uncertainties: uncertainties in the solubility in log(S)/L
    """

    solubilities: Annotated[list[float], AfterValidator(round_list(6))]
    uncertainties: Annotated[list[float | None], AfterValidator(round_list(6))]

    @model_validator(mode="after")
    def check_size(self) -> Self:
        """Confirm that the shapes are the same"""
        if len(self.solubilities) != len(self.uncertainties):
            raise ValueError("Solubilities and uncertainties must have the same length.")

        return self


class SolubilityWorkflow(SMILESWorkflow):
    """
    Solubility prediction workflow.

    Inputs:
    :param initial_smiles: SMILES string of the molecule
    :param solubility_method: model used for solubility prediction
    :param solvents: list of solvent SMILES strings
    :param temperatures: temperatures in K

    Results:
    :param solubilities: {solvent: SolubilityResult}
    """

    initial_smiles: str
    solubility_method: SolubilityMethod = SolubilityMethod.FASTSOLV
    solvents: list[str] = ["O"]
    temperatures: list[float] = [298.15]

    solubilities: dict[str, SolubilityResult] = {}

    @model_validator(mode="after")
    def check_size(self) -> Self:
        """Check that the sizes of the lists are consistent."""
        for solvent, result in self.solubilities.items():
            if len(result.solubilities) != len(self.temperatures):
                raise ValueError(f"Solubilities for {solvent} must have the same length as temperatures.")

            if solvent not in self.solvents:
                raise ValueError(f"Solvent {solvent} not in initial list of solvents.")

        return self

    @model_validator(mode="after")
    def check_solvent_temperature(self) -> Self:
        """Check that models with limited domain of applicability are predicting within correct domain."""
        match self.solubility_method:
            case SolubilityMethod.KINGFISHER | SolubilityMethod.ESOL:
                if (len(self.solvents) > 1) or (self.solvents[0] != "O"):
                    raise ValueError(f"Method `{self.solubility_method}` can only predict aqueous solubility, so `solvents` must be [`O`] only.")
                if (len(self.temperatures) > 1) or (abs(self.temperatures[0] - 298.15) > 0.1):
                    raise ValueError(
                        f"Method `{self.solubility_method}` can only predict solubility at room temperature, so `temperatures` must be [298.15] only."
                    )
            case _:
                pass

        return self
