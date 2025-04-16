"""Basic calculation workflow."""

from typing import Any

from ..types import UUID
from .workflow import FastaWorkflow


class StructurePredictionWorkflow(FastaWorkflow):
    """
    A workflow for predicting structures. Especially protein structures.

    Inherited:
    :param fasta_string: fasta file string to run stucture prediction on

    New:
    :param use_msa_server: whether to use the MSA server
    :param use_templates_server: whether to use the templates server
    :param predicted_structure_uuid: UUID of the predicted structure
    """

    use_msa_server: bool = False
    use_templates_server: bool = False
    predicted_structure_uuid: UUID | None = None
    prediciton_scores: dict[str, Any] | None = None
