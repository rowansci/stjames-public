"""Ion mobility workflow."""

from ..types import UUID
from .workflow import MoleculeWorkflow


class IonMobilityWorkflow(MoleculeWorkflow):
    """
    Workflow for calculating hydrogen bond basicity.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param do_csearch: whether to perform a conformational search
    :param do_optimization: whether to perform an optimization

    Results:
    :param conformer_ccs: the collision cross section (Å**2) per conformer
    :param conformer_ccs_stdev: the uncertainty in the same
    :param conformer_weights: the Boltzmann weights at RT
    :param average_ccs: the Boltzmann-weighted CCS for the ensemble
    :param average_ccs_stdev: the uncertainty in the same
    """

    do_csearch: bool = True
    do_optimization: bool = True
    conformers: list[UUID] = []

    conformer_ccs: list[float] = []
    conformer_ccs_stdev: list[float] = []
    boltzmann_weights: list[float] = []

    average_ccs: float | None = None
    average_ccs_stdev: float | None = None
