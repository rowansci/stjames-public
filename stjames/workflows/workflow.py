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


class FASTAWorkflow(Workflow):
    """
    Base class for Workflows that operate on protein sequences and SMILES.

    :param initial_protein_sequences: protein sequences of interest
    :param initial_smiles_list: SMILES strings of interest
    """

    initial_protein_sequences: list[str]
    initial_smiles_list: list[str] | None = None
    ligand_binding_affinity_index: int | None = None

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.initial_protein_sequences} {self.initial_smiles_list}>"


class SMILESWorkflow(Workflow):
    """
    Base class for Workflows that operate on a single SMILES string.

    :param initial_smiles: SMILES string of interest
    """

    initial_smiles: str

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.initial_smiles}>"


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
