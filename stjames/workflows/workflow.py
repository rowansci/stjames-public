"""Base classes for workflows."""

from pydantic import field_validator, model_validator

from ..base import Base
from ..message import Message
from ..mode import Mode
from ..molecule import Molecule
from ..types import UUID
from .dna import DNASequence
from .protein import ProteinSequence
from .rna import RNASequence


class Workflow(Base):
    """
    Base class for Workflows.

    :param messages: messages to display
    """

    messages: list[Message] = []

    def __str__(self) -> str:
        return repr(self)


class FASTAWorkflow(Workflow):
    """
    Base class for Workflows that operate on biological sequences and SMILES.

    :param initial_protein_sequences: protein sequences to evaluate, either plain sequence strings or ProteinSequence objects with metadata (optional)
    :param initial_dna_sequences: DNA sequences to evaluate, either plain sequence strings or DNASequence objects with metadata (optional)
    :param initial_rna_sequences: RNA sequences to evaluate, either plain sequence strings or RNASequence objects with metadata (optional)
    :param initial_smiles_list: SMILES strings of interest
    :param ligand_binding_affinity_index: optional index selecting which ligand affinity to evaluate
    :raises ValueError: if none of the sequence lists are provided
    """

    initial_protein_sequences: list[ProteinSequence | str] | None = None
    initial_dna_sequences: list[DNASequence | str] | None = None
    initial_rna_sequences: list[RNASequence | str] | None = None
    initial_smiles_list: list[str] | None = None
    ligand_binding_affinity_index: int | None = None

    @model_validator(mode="after")
    @classmethod
    def _require_sequence(cls, model: "FASTAWorkflow") -> "FASTAWorkflow":
        """Ensure at least one biological sequence list is populated."""

        sequence_lists = (
            model.initial_protein_sequences,
            model.initial_dna_sequences,
            model.initial_rna_sequences,
        )

        if not any(seq for seq in sequence_lists if seq):
            raise ValueError(
                "Provide at least one of `initial_protein_sequences`, `initial_dna_sequences`, or `initial_rna_sequences`.",
            )

        return model


class SMILESWorkflow(Workflow):
    """
    Base class for Workflows that operate on a single SMILES string.

    :param initial_smiles: SMILES string of interest
    """

    initial_smiles: str

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.initial_smiles}>"


class BatchSMILESWorkflow(Workflow):
    """
    Base class for Workflows that operate on a list of SMILES strings.

    :param smiles_list: SMILES strings of interest
    """

    smiles_list: list[str]


class MoleculeWorkflow(Workflow):
    """
    Base class for Workflows that operate on a single molecule.

    :param initial_molecule: Molecule of interest
    :param mode: Mode to use
    """

    initial_molecule: Molecule
    mode: Mode = Mode.AUTO

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.mode.name}>"

    @field_validator("mode")
    @classmethod
    def set_mode_auto(cls, mode: Mode) -> Mode:
        """Set the mode to RAPID if AUTO is selected."""
        if mode == Mode.AUTO:
            return Mode.RAPID

        return mode


class ProteinSequenceWorkflow(Workflow):
    """
    Base class for Workflows that operate on protein sequences.

    :param protein_sequences: protein sequences of interest
    """

    protein_sequences: list[str]


class DBCalculation(Base):
    """Encodes a calculation that's in the database. This isn't terribly useful by itself."""

    uuid: UUID
