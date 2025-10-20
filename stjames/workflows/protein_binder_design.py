"""Protein-binder-design workflow."""

from typing import TypeAlias

from ..base import Base, LowercaseStrEnum
from ..types import UUID
from .protein_cofolding import AffinityScore, CofoldingScores, Token
from .workflow import FASTAWorkflow

ProteinUUID: TypeAlias = UUID


class ProteinBinderDesignResult(Base):
    """
    The output; a designed binder.

    :param sequence: the sequence
    :param isolated_structure: the PDB of the structure on its own
    :param bound_structure: the PDB of the structure bound to the target
    :param affinity_score: the predicted affinity
    :param scores: the co-folding scores for the generated structure
    """

    sequence: str
    isolated_structure: ProteinUUID | None = None
    bound_structure: ProteinUUID | None = None
    affinity_score: AffinityScore | None = None
    scores: CofoldingScores | None = None


class BondConstraint(Base):
    """
    Specifies bonds between two chains.
    """

    token1: Token
    token2: Token


class BoltzGenProtocol(LowercaseStrEnum):
    """
    The predefined protocol used for generation + filtering.
    """

    PROTEIN_ANYTHING = "protein-anything"
    PEPTIDE_ANYTHING = "peptide-anything"
    PROTEIN_SMALL_MOLECULE = "protein-small_molecule"
    NANOBODY_ANYTHING = "nanobody-anything"


class BoltzGenSettings(Base):
    """
    The settings for running BoltzGen.

    :param num_designs: how many designs to generate
    :param protocol: which protocol to use
    :param binding_residue: a dict mapping the chain ID to which residues should bind.
        the string follows the BoltzGen format of specifying ranges of residue indices (refer to their documentation).
        examples include "5..7,13" or "5..15,50..".
    """

    protocol: BoltzGenProtocol = BoltzGenProtocol.PROTEIN_ANYTHING
    num_designs: int = 100
    binding_residues: dict[int, str] = {}


class ProteinBinderDesignWorkflow(FASTAWorkflow):
    """
    A workflow for generating proteins or peptides that bind to something.

    Inherited:
    :param initial_protein_sequences: protein sequences of interest
    :param initial_smiles_list: SMILES strings of interest

    New:
    :param input_protein_uuid: the input protein structure, if requested
    :param binder_design_settings: the settings for the protein generation method employed
    :param bond_constraints: BoltzGen bond constraints
    :param generated_binders: the output structures
    """

    input_protein_uuid: ProteinUUID | None = None

    binder_design_settings: BoltzGenSettings = BoltzGenSettings()
    bond_constraints: list[BondConstraint] = []

    generated_binders: list[ProteinBinderDesignResult] = []
