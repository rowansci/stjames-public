"""DNA-related workflow data models."""

from ..base import Base


class DNASequence(Base):
    """
    DNA sequence metadata.

    :param sequence: nucleotide sequence string
    """

    sequence: str
