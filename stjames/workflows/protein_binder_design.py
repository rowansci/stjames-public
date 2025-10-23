"""Protein-binder-design workflow."""

from typing import TypeAlias

from ..base import Base, LowercaseStrEnum
from ..types import UUID
from .protein_cofolding import AffinityScore, Token
from .workflow import FASTAWorkflow

ProteinUUID: TypeAlias = UUID


class ProteinBinderScores(Base):
    """
    Compact, interpretable metrics for a designed binder.
    ↑ higher is better, ↓ lower is better

    :param quality_score: aggregate model quality (↑)
    :param num_filters_passed: number of QC/heuristic filters passed (↑)
    :param iptm: inter-chain pTM confidence, 0–1 (↑)
    :param design_ptm: design pTM confidence, 0–1 (↑)
    :param min_interaction_pae: minimum interface PAE in Å (↓)
    :param bb_rmsd: backbone RMSD in Å (↓)
    :param delta_sasa_refolded: ΔSASA of interface after refolding, Å² (↑ typically indicates better burial)
    :param plip_hbonds_refolded: count of hydrogen bonds at the interface (↑)
    :param plip_saltbridge_refolded: count of salt bridges at the interface (↑)
    :param liability_score: composite liabilities score (↓)
    :param liability_high_severity_violations: count of high-severity liabilities (↓)
    :param liability_num_violations: total liability count (↓)
    :param helix: fraction helical content, 0–1
    :param sheet: fraction β-sheet content, 0–1
    :param loop: fraction loop/coil content, 0–1
    :param design_largest_hydrophobic_patch_refolded: largest hydrophobic patch area after refolding, Å²
    :param design_hydrophobicity: overall design hydrophobicity score (unitless)
    :param num_tokens: sequence length / token count
    """

    quality_score: float | None = None
    num_filters_passed: int | None = None

    iptm: float | None = None
    design_ptm: float | None = None
    min_interaction_pae: float | None = None

    bb_rmsd: float | None = None
    delta_sasa_refolded: float | None = None

    plip_hbonds_refolded: int | None = None
    plip_saltbridge_refolded: int | None = None

    liability_score: float | None = None
    liability_high_severity_violations: int | None = None
    liability_num_violations: int | None = None

    helix: float | None = None
    sheet: float | None = None
    loop: float | None = None
    design_largest_hydrophobic_patch_refolded: float | None = None
    design_hydrophobicity: float | None = None
    num_tokens: int | None = None


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
    scores: ProteinBinderScores | None = None


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
    budget: int = 20
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

    input_protein_uuids: list[ProteinUUID] = []

    binder_design_settings: BoltzGenSettings = BoltzGenSettings()
    bond_constraints: list[BondConstraint] = []

    generated_binders: list[ProteinBinderDesignResult] = []
