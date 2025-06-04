"""Hydrogen-bond-basicity workflow."""

from ..base import Base
from ..types import UUID
from .workflow import MoleculeWorkflow


class HydrogenBondAcceptorSite(Base):
    """
    A hydrogen-bond-acceptor site.

    :param atom_idx: index of the atom
    :param pkbhx: Hydrogen-bond basicity
    :param position: position of the atom
    :param name: name of the atom
    """

    atom_idx: int  # zero-indexed
    pkbhx: float
    position: tuple[float, float, float]
    name: str | None = None


class HydrogenBondDonorSite(Base):
    """
    A hydrogen-bond-donor site.

    :param atom_idx: index of the atom
    :param pk_alpha: Hydrogen-bond acidity
    :param position: position of the atom
    """

    atom_idx: int  # zero-indexed
    pk_alpha: float
    position: tuple[float, float, float]


class HydrogenBondBasicityWorkflow(MoleculeWorkflow):
    """
    Workflow for calculating hydrogen-bond basicity and acidity.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param do_csearch: whether to perform a conformational search
    :param do_optimization: whether to perform an optimization

    Results:
    :param optimization: UUID of optimization
    :param hba_sites: hydrogen-bond-acceptor sites
    :param hbd_sites: hydrogen-bond-donor sites
    """

    do_csearch: bool = True
    do_optimization: bool = True

    optimization: UUID | None = None
    hba_sites: list[HydrogenBondAcceptorSite] = []  # noqa: RUF012
    hbd_sites: list[HydrogenBondDonorSite] = []  # noqa: RUF012
