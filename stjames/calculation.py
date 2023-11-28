from typing import Optional

from .molecule import Molecule
from .base import Base, LowercaseStrEnum
from .settings import Settings
from .status import Status


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

    engine: Optional[str] = "kestrel"

    # not to be changed by end users, diff. versions will have diff. defaults
    json_format: str = StJamesVersion.V0
