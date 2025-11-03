"""RNA-related data models."""

from .base import Base


class RNASequence(Base):
    """
    RNA sequence data.

    :param sequence: nucleotide string
    """

    sequence: str
