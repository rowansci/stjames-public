"""Protein Cofolding Workflow."""

from pydantic import BaseModel

from ..base import LowercaseStrEnum
from ..types import UUID
from .workflow import FASTAWorkflow


class CofoldingModel(LowercaseStrEnum):
    """Cofolding model to be used for prediction."""

    CHAI_1R = "chai_1r"
    BOLTZ_1 = "boltz_1"
    BOLTZ_2 = "boltz_2"


class CofoldingScores(BaseModel):
    confidence_score: float
    ptm: float  # predicted template modeling score
    iptm: float  # interface predicted template modeling score
    avg_lddt: float


class AffinityScore(BaseModel):
    pred_value: float
    probability_binary: float
    pred_value1: float
    probability_binary1: float
    pred_value2: float
    probability_binary2: float


class ProteinCofoldingWorkflow(FASTAWorkflow):
    """
    A workflow for predicting structures. Especially protein structures.

    Inherited:
    :param initial_protein_sequences: protein sequences of interest
    :param initial_smiles_list: SMILES strings of interest

    New:
    :param use_msa_server: whether to use the MSA server
    :param use_templates_server: whether to use the templates server
    :param predicted_structure_uuid: UUID of the predicted structure
    """

    use_msa_server: bool = False
    use_templates_server: bool = False
    use_potentials: bool = False
    predicted_structure_uuid: UUID | None = None
    scores: CofoldingScores | None = None
    model: CofoldingModel = CofoldingModel.BOLTZ_2
    affinity_score: AffinityScore | None = None
    lddt: list[float] | None = None
