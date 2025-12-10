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
    :param caco_2_P_app: Caco-2 apparent permeability P_app (x10^-6 cm/s)
    """

    initial_smiles: str
    membrane_permeability_method: MembranePermeabilityMethod = MembranePermeabilityMethod.CHEMPROP_OHLSSON2025

    caco_2_P_app: Annotated[float, AfterValidator(round_float(3))] | None = None
