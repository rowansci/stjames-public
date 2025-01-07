from ..base import Base
from ..types import UUID
from .workflow import Workflow


class HydrogenBondAcceptorSite(Base):
    atom_idx: int  # zero-indexed
    pkbhx: float
    position: tuple[float, float, float]
    name: str | None = None


class HydrogenBondBasicityWorkflow(Workflow):
    do_csearch: bool = True
    do_optimization: bool = True

    # UUID of optimization
    optimization: UUID | None = None

    # hydrogen-bond-acceptor sites
    hba_sites: list[HydrogenBondAcceptorSite] = []  # noqa: RUF012
