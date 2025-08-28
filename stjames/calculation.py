from typing import Optional, Self

from pydantic import model_validator

from .base import Base, LowercaseStrEnum, UniqueList
from .message import Message
from .molecule import Molecule
from .settings import Settings
from .status import Status
from .task import Task
from .types import UUID


class StJamesVersion(LowercaseStrEnum):
    """
    Flag for which version of ``stjames`` this calculation was saved in.
    This will become meaningful later, when we have (probably) multiple versions.
    """

    V0 = "09072023"


class Calculation(Base):
    molecules: list[Molecule]

    tasks: UniqueList[Task] = []
    settings: Settings = Settings()

    status: Status = Status.QUEUED

    name: Optional[str] = None
    elapsed: Optional[float] = None
    logfile: Optional[str] = None
    messages: list[Message] = []

    # DEPRECATED - moving into settings
    engine: Optional[str] = "peregrine"

    uuids: list[UUID | None] | None = None

    # not to be changed by end users, diff. versions will have diff. defaults
    json_format: str = StJamesVersion.V0

    @model_validator(mode="after")
    def populate_tasks(self) -> Self:
        """Set the tasks from the settings, so that we don't have to migrate old entries."""
        if len(self.tasks) == 0:
            self.tasks = self.settings.tasks
        return self
