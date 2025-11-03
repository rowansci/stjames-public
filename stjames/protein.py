"""Protein-related data models."""

from .base import Base


class ProteinSequence(Base):
    """
    Protein sequence metadata including cyclic flag.

    :param sequence: amino-acid sequence string
    :param cyclic: whether this sequence forms a cyclic peptide (defaults to False)
    """

    sequence: str
    cyclic: bool = False
