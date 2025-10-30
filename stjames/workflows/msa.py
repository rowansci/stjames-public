from ..base import LowercaseStrEnum
from .workflow import ProteinSequenceWorkflow


class MSAFormat(LowercaseStrEnum):
    """Format of the MSA."""

    COLABFOLD_DEFAULT = "colabfold_default"
    AF3_JSON = "af3_json"


class MSAWorkflow(ProteinSequenceWorkflow):
    """
    Workflow for generating a MSA from protein sequences.

    Inherited:
    :param initial_protein_sequences: protein sequences of interest

    New:
    :param format: the format of the MSA return files

    Results:
    :param a3m_file: A3M file string
    :param m8_file: M8 file string
    :param af3_json_file: AF3 JSON file string
    """

    format: MSAFormat = MSAFormat.COLABFOLD_DEFAULT
    a3m_file: str | None = None
    m8_file: str | None = None
    af3_json_file: str | None = None
