"""pKa workflow."""

from typing import Annotated, Optional, Self

from pydantic import AfterValidator, model_validator

from ..base import Base, round_float
from ..types import round_list
from .workflow import SMILESWorkflow


class MacropKaMicrostate(Base):
    """
    A microstate for pKa calculations.

    :param smiles: SMILES string for this conformer
    :param energy: free energy of this conformer
    :param charge: the total charge
    """

    smiles: str
    energy: Annotated[float, AfterValidator(round_float(3))]  # free energy
    charge: int


class MacropKaValue(Base):
    """
    Represents a change in pKa.

    :param initial_charge: the charge of the initial state
    :param final_charge: the charge of the final state
    :param pKa: the pKa for the transition
    """

    initial_charge: int
    final_charge: int
    pKa: Annotated[float, AfterValidator(round_float(3))]


class MacropKaWorkflow(SMILESWorkflow):
    """
    Workflow for calculating pKa.

    Inherited:
    :param initial_smiles:

    New:
    :param temperature: the temperature, in K
    :param min_pH: for precomputed microstate weights
    :param max_pH: for precomputed microstate weights

    Results:
    :param microstates: microstates
    :param pKa_values: macroscopic pKa values
    :param microstate_weights_by_pH: precompute the % of different microstates
    """

    min_pH: Annotated[float, AfterValidator(round_float(3))] = 0.0
    max_pH: Annotated[float, AfterValidator(round_float(3))] = 14.0

    max_charge: int = 2
    min_charge: int = -2

    microstates: list[MacropKaMicrostate] = []
    pKa_values: list[MacropKaValue] = []
    microstate_weights_by_pH: list[
        tuple[
            Annotated[float, AfterValidator(round_float(3))],
            Annotated[list[float], AfterValidator(round_list(6))],
        ]
    ] = []

    isoelectric_point: Annotated[Optional[float], AfterValidator(round_float(3))] = None

    logD_by_pH: list[
        tuple[
            Annotated[float, AfterValidator(round_float(3))],
            Annotated[float, AfterValidator(round_float(3))],
        ]
    ] = []

    aqueous_solubility_by_pH: list[
        tuple[
            Annotated[float, AfterValidator(round_float(3))],
            Annotated[float, AfterValidator(round_float(3))],
        ]
    ] = []

    @model_validator(mode="after")
    def check_weights(self) -> Self:
        for _, weights in self.microstate_weights_by_pH:
            if len(weights) != len(self.microstates):
                raise ValueError("Length of microstate weights doesn't match!")

        return self

    @model_validator(mode="after")
    def check_minmax_charges(self) -> Self:
        if self.min_charge >= self.max_charge:
            raise ValueError("Incoherent min/max charge specification")

        return self
