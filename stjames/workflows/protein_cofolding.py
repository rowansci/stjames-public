"""Protein Cofolding Workflow."""

from pydantic import BaseModel

from ..types import UUID
from .workflow import FASTAWorkflow


class CofoldingScores(BaseModel):
    confidence_score: float
    ptm: float  # predicted template modeling score
    iptm: float  # interface predicted template modeling score


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
    predicted_structure_uuid: UUID | None = None
    scores: CofoldingScores | None = None
