"""Solubility prediction workflow."""

from typing import Annotated, Self

from pydantic import AfterValidator, BaseModel, model_validator

from ..types import round_list


class SolubilityResult(BaseModel):
    """
    Solubility result.

    :param solubilities: solubilities in log(S)/L
    :param uncertainties: uncertainties in the solubility in log(S)/L
    """

    solubilities: Annotated[list[float], AfterValidator(round_list(6))]
    uncertainties: Annotated[list[float], AfterValidator(round_list(6))]

    @model_validator(mode="after")
    def check_size(self) -> Self:
        """Confirm that the shapes are the same"""
        if len(self.solubilities) != len(self.uncertainties):
            raise ValueError("Solubilities and uncertainties must have the same length.")

        return self


class SolubilityWorkflow(BaseModel):
    """
    Solubility prediction workflow.

    Inputs:
    :param initial_smiles: SMILES string of the molecule
    :param solvents: list of solvent SMILES strings
    :param temperatures: list of temperatures in K

    Results:
    :param solubilities: {solvent: SolubilityResult}
    """

    initial_smiles: str
    solvents: list[str]
    temperatures: list[float]

    solubilities: dict[str, SolubilityResult] = {}

    @model_validator(mode="after")
    def check_size(self) -> Self:
        """Check that the sizes of the lists are consistent."""
        for solvent, result in self.solubilities.items():
            if len(result.solubilities) != len(self.temperatures):
                raise ValueError(f"Solubilities for {solvent} must have the same length as temperatures.")

            if solvent not in self.solvents:
                raise ValueError(f"Solvent {solvent} not in initial list of solvents")

        return self
