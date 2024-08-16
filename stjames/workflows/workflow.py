from pydantic import BaseModel, ConfigDict

from ..base import Base
from ..message import Message
from ..mode import Mode
from ..molecule import Molecule
from ..types import UUID


class Workflow(Base):
    """All workflows should have these properties."""

    initial_molecule: Molecule
    messages: list[Message] = []


class DBCalculation(Base):
    """Encodes a calculation that's in the database. This isn't terribly useful by itself."""

    uuid: UUID


class WorkflowInput(BaseModel):
    """
    Input for a workflow.

    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow
    """

    model_config = ConfigDict(extra="forbid")

    initial_molecule: Molecule
    mode: Mode

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f"<{type(self).__name__} {self.mode.name}>"


class WorkflowResults(BaseModel):
    """Results of a workflow."""

    model_config = ConfigDict(extra="forbid", frozen=True)
