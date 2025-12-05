from typing import Annotated

from pydantic import AfterValidator

from ..base import LowercaseStrEnum, round_float
from .workflow import SMILESWorkflow


class MembranePermeabilityMethod(LowercaseStrEnum):
    CHEMPROP_OHLSSON2025 = "chemprop_ohlsson2025"


class MembranePermeabilityWorkflow(SMILESWorkflow):
    """
    Membrane permeability prediction workflow.

    Inputs:
    :param initial_smiles: SMILES string of the molecule
    :param membrane_permeability_method: model used to predict membrane permeability

    Results:
    :param membrane_permeability: the log(P_app) in cm/s of the compound across a membrane
    """

    initial_smiles: str
    membrane_permeability_method: MembranePermeabilityMethod = MembranePermeabilityMethod.CHEMPROP_OHLSSON2025

    membrane_permeability: Annotated[float, AfterValidator(round_float(3))]
