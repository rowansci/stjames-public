from typing import Optional

from ..base import Base
from ..method import Method
from ..solvent import Solvent
from .workflow import Workflow


class ConformerSettings(Base):
    csearch_program: str = "rdkit"

    num_confs_considered: int = 100
    num_confs_taken: int = 50

    final_method: Method = Method.AIMNET2_WB97MD3
    solvent: Optional[Solvent] = Solvent.WATER
    max_energy: float = 5


class ETKDGConformerSettings(ConformerSettings):
    num_initial_confs: int = 100
    max_mmff_energy: float = 10
    num_confs_considered: float = 10
    num_confs_taken: float = 3
    rmsd_cutoff: float = 0.25


class CrestConformerSettings(ConformerSettings):
    flags: str = "--quick --ewin 10"
    gfn: str = "ff"


class Conformer(Base):
    energy: float
    weight: Optional[float] = None

    # uuid, optionally
    uuid: Optional[str] = None


class ConformerWorkflow(Workflow):
    settings: ConformerSettings
    conformers: list[Conformer] = []
