from ..base import Base
from ..molecule import Molecule


class Workflow(Base):
    """All workflows should have these properties."""

    initial_molecule: Molecule


class DBCalculation(Base):
    """Encodes a calculation that's in the database. This isn't terribly useful by itself."""

    uuid: str
