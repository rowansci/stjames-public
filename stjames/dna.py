"""DNA-related data models."""

from .base import Base


class DNASequence(Base):
    """
    DNA sequence data.

    :param sequence: nucleotide string
    """

    sequence: str
