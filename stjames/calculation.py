from typing import Optional

from .base import Base, LowercaseStrEnum
from .message import Message
from .molecule import Molecule
from .settings import Settings
from .status import Status
from .types import UUID


class StJamesVersion(LowercaseStrEnum):
    """
    Flag for which version of ``stjames`` this calculation was saved in.
    This will become meaningful later, when we have (probably) multiple versions.
    """

    V0 = "09072023"


class Calculation(Base):
    molecules: list[Molecule]

    settings: Settings = Settings()

    status: Status = Status.QUEUED

    name: Optional[str] = None
    elapsed: Optional[float] = None
    logfile: Optional[str] = None
    messages: list[Message] = []

    engine: Optional[str] = "peregrine"
    uuids: list[UUID | None] | None = None

    # not to be changed by end users, diff. versions will have diff. defaults
    json_format: str = StJamesVersion.V0
