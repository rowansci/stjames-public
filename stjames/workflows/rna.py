"""RNA-related workflow data models."""

from ..base import Base


class RNASequence(Base):
    """
    RNA sequence metadata.

    :param sequence: nucleotide string
    """

    sequence: str
