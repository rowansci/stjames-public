from ..base import LowercaseStrEnum
from .workflow import ProteinSequenceWorkflow


class MSAFormat(LowercaseStrEnum):
    """Format of the MSA."""

    COLABFOLD = "colabfold"
    CHAI = "chai"
    BOLTZ = "boltz"


class MSAWorkflow(ProteinSequenceWorkflow):
    """
    Workflow for generating a MSA from protein sequences.

    Inherited:
    :param initial_protein_sequences: protein sequences of interest

    New:
    :param output_formats: the formats of the MSA return files
    """

    output_formats: list[MSAFormat] = [MSAFormat.COLABFOLD]
