"""DNA-related tdata models."""

from .base import Base


class DNASequence(Base):
    """
    DNA sequence metadata.

    :param sequence: nucleotide string
    """

    sequence: str
