from ..base import Base
from ..mode import Mode
from ..molecule import Molecule


class Workflow(Base):
    """All workflows should have these properties."""

    initial_molecule: Molecule
    mode: Mode


class DBCalculation(Base):
    """Encodes a calculation that's in the database. This isn't terribly useful by itself."""

    uuid: str
