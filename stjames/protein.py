"""Protein-related data models."""

from .base import Base


class ProteinSequence(Base):
    """
    Protein sequence data.

    :param sequence: amino-acid sequence string
    :param cyclic: whether this sequence forms a cyclic peptide
    """

    sequence: str
    cyclic: bool = False
