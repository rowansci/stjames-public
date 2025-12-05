from typing import Annotated

from pydantic import AfterValidator

from ..base import LowercaseStrEnum, round_float
from .workflow import SMILESWorkflow


class CellMembranePermeabilityMethod(LowercaseStrEnum):
    CHEMPROP_OHLSSON2025 = "chemprop_ohlsson2025"


class CellMembranePermeabilityWorkflow(SMILESWorkflow):
    """
    Cell membrane permeability prediction workflow.

    Inputs:
    :param initial_smiles: SMILES string of the molecule
    :param cell_membrane_permeability_method: model used to predict cell membrane permeability

    Results:
    :param cell_membrane_permeability: the log(P_app) in cm/s of the compound across a cell membrane
    """

    initial_smiles: str
    cell_membrane_permeability_method: CellMembranePermeabilityMethod = CellMembranePermeabilityMethod.CHEMPROP_OHLSSON2025

    cell_membrane_permeability: Annotated[float, AfterValidator(round_float(3))]
