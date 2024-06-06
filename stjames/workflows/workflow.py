from ..base import Base
from ..mode import Mode
from ..molecule import Molecule


class Workflow(Base):
    initial_molecule: Molecule
    mode: Mode
