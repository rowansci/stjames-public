from typing import Annotated, Self

from pydantic import AfterValidator, model_validator

from ..base import LowercaseStrEnum, round_float
from ..molecule import Molecule
from .workflow import MoleculeWorkflow, SMILESWorkflow


class MembranePermeabilityMethod(LowercaseStrEnum):
    CHEMPROP_OHLSSON2025 = "chemprop_ohlsson2025"
    PYPERMM = "pypermm"


class MembranePermeabilityWorkflow(SMILESWorkflow, MoleculeWorkflow):
    """
    Membrane permeability prediction workflow.

    Inherited:
    :param initial_smiles: SMILES string of the molecule
    :param initial_molecule: SMILES string of the molecule

    Inputs:
    :param membrane_permeability_method: model used to predict membrane permeability

    Results:
    :param caco_2_P_app: base-10 logarithm of the Caco-2 apparent permeability P_app (in cm/s)
    :param caco_2_logP: the Caco-2 intrinsic permeability coefficient logP
    :param blm_logP: the bilayer-lipid-membrane intrinsic permeability coefficient logP
    :param pampa_logP: the PAMPA intrinsic permeability coefficient logP
    :param plasma_logP: the plasma-membrane intrinsic permeability coefficient logP
    :param bbb_logP: the blood–brain-barrier intrinsic permeability coefficient logP
    :param energy_profile: how the energy of the compound is predicted to change as it passes through a membrane.
        the first value is the membrane profile (Å), the second is the energy (kcal/mol).
    """

    initial_smiles: str = ""
    initial_molecule: Molecule | None = None  # type: ignore [assignment]

    membrane_permeability_method: MembranePermeabilityMethod = MembranePermeabilityMethod.CHEMPROP_OHLSSON2025

    caco_2_P_app: Annotated[float, AfterValidator(round_float(3))] | None = None

    caco_2_logP: Annotated[float, AfterValidator(round_float(3))] | None = None
    blm_logP: Annotated[float, AfterValidator(round_float(3))] | None = None
    pampa_logP: Annotated[float, AfterValidator(round_float(3))] | None = None
    plasma_logP: Annotated[float, AfterValidator(round_float(3))] | None = None
    bbb_logP: Annotated[float, AfterValidator(round_float(3))] | None = None

    energy_profile: list[
        tuple[
            Annotated[float, AfterValidator(round_float(3))],
            Annotated[float, AfterValidator(round_float(3))],
        ]
    ] = []

    @model_validator(mode="after")
    def validate_mol_input(self) -> Self:
        """Ensure that only one of initial_molecule or initial_smiles is set."""

        if not (bool(self.initial_smiles) ^ bool(self.initial_molecule)):
            raise ValueError("Can only set one of initial_molecule and initial_smiles")

        if self.initial_molecule:
            assert self.membrane_permeability_method == MembranePermeabilityMethod.PYPERMM
        elif self.initial_smiles:
            assert self.membrane_permeability_method == MembranePermeabilityMethod.CHEMPROP_OHLSSON2025

        return self
