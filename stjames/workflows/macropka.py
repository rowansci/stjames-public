"""pKa workflow."""

from typing import Annotated, Optional, Self

from pydantic import AfterValidator, model_validator

from ..base import Base, round_float, round_optional_float
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
    :param min_pH: for precomputed microstate weights, logD, etc
    :param max_pH: for precomputed microstate weights, logD, etc
    :param max_charge: max charge to consider for microstates
    :param min_charge: min charge to consider for microstates
    :param compute_aqueous_solubility: whether or not to compute aqueous solubility
    :param compute_solvation_energy: whether to run a csearch + compute the solvation energy (for Kpuu)

    Results:
    :param microstates: microstates
    :param pKa_values: macroscopic pKa values
    :param isoelectric_point: the isoelectric point (in pH units)
    :param solvation_energy: the solvation energy, in kcal/mol
    :param microstate_weights_by_pH: the % of different microstates by pH
    :param logD_by_pH: the distribution constant (water/octanol) by pH
    :param aqueous_solubility_by_pH: the log(S)/L of the compound in water, by pH
    :param kpuu_probability: the probability that Kpuu >= 0.3, the Schrodinger-determined threshold
    """

    min_pH: Annotated[float, AfterValidator(round_float(3))] = 0.0
    max_pH: Annotated[float, AfterValidator(round_float(3))] = 14.0
    max_charge: int = 2
    min_charge: int = -2
    compute_aqueous_solubility: bool = False
    compute_solvation_energy: bool = True

    microstates: list[MacropKaMicrostate] = []
    pKa_values: list[MacropKaValue] = []
    isoelectric_point: Annotated[Optional[float], AfterValidator(round_optional_float(3))] = None
    solvation_energy: Annotated[Optional[float], AfterValidator(round_optional_float(3))] = None
    kpuu_probability: Annotated[Optional[float], AfterValidator(round_optional_float(3))] = None

    microstate_weights_by_pH: list[
        tuple[
            Annotated[float, AfterValidator(round_float(3))],
            Annotated[list[float], AfterValidator(round_list(6))],
        ]
    ] = []

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
