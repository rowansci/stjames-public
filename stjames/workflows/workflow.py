"""Base classes for workflows."""

from pydantic import field_validator

from ..base import Base
from ..message import Message
from ..mode import Mode
from ..molecule import Molecule
from ..types import UUID


class Workflow(Base):
    """
    Base class for Workflows.

    :param messages: messages to display
    """

    messages: list[Message] = []

    def __str__(self) -> str:
        return repr(self)


class ProteinSequence(Base):
    """
    Protein sequence metadata including cyclic flag.

    :param sequence: amino-acid sequence string
    :param cyclic: whether this sequence forms a cyclic peptide (defaults to False)
    """

    sequence: str
    cyclic: bool = False


class FASTAWorkflow(Workflow):
    """
    Base class for Workflows that operate on protein sequences and SMILES.

    :param initial_protein_sequences: proteins to evaluate, either plain sequence strings or ProteinSequence objects with cyclic flags
    :param initial_smiles_list: SMILES strings of interest
    :param ligand_binding_affinity_index: optional index selecting which ligand affinity to evaluate
    """

    initial_protein_sequences: list[ProteinSequence] | list[str]
    initial_smiles_list: list[str] | None = None
    ligand_binding_affinity_index: int | None = None


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


class DBCalculation(Base):
    """Encodes a calculation that's in the database. This isn't terribly useful by itself."""

    uuid: UUID
