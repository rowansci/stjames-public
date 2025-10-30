from ..base import LowercaseStrEnum
from .workflow import ProteinSequenceWorkflow


class MsaFormat(LowercaseStrEnum):
    """Format of the MSA."""

    COLABFOLD_DEFAULT = "colabfold_default"
    AF3_JSON = "af3_json"


class MsaWorkflow(ProteinSequenceWorkflow):
    """
    Workflow for generating a MSA from protein sequences.

    Inherited:
    :param initial_protein_sequences: protein sequences of interest

    New:
    :param format: the format of the MSA return files

    """

    format: MsaFormat = MsaFormat.COLABFOLD_DEFAULT
