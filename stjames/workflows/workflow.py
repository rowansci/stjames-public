from pydantic import field_validator

from ..base import Base
from ..message import Message
from ..mode import Mode
from ..molecule import Molecule
from ..types import UUID


class Workflow(Base):
    """
    Base class for Workflows.

    :param initial_molecule: Molecule of interest
    :param mode: Mode to use
    :param messages: messages to display
    """

    initial_molecule: Molecule
    mode: Mode = Mode.AUTO
    messages: list[Message] = []

    def __str__(self) -> str:
        return repr(self)

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
