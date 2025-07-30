"""Ion mobility workflow."""

from typing import Annotated

from pydantic import AfterValidator

from ..base import round_optional_float
from ..types import UUID, round_list
from .workflow import MoleculeWorkflow


class IonMobilityWorkflow(MoleculeWorkflow):
    """
    Workflow for calculating hydrogen bond basicity.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param protonate: automatically protonate the molecule
    :param temperature: the temperature, in Kelvin
    :param do_csearch: whether to perform a conformational search
    :param do_optimization: whether to perform an optimization

    Results:
    :param conformers: the UUIDs of the conformers
    :param conformer_ccs: the collision cross section (Å**2) per conformer
    :param conformer_ccs_stdev: the uncertainty in the same
    :param conformer_weights: the Boltzmann weights at RT
    :param average_ccs: the Boltzmann-weighted CCS for the ensemble
    :param average_ccs_stdev: the uncertainty in the same
    """

    protonate: bool = False
    temperature: float = 300
    do_csearch: bool = True
    do_optimization: bool = True

    conformers: list[UUID] = []

    conformer_ccs: Annotated[list[float], AfterValidator(round_list(3))] = []
    conformer_ccs_stdev: Annotated[list[float], AfterValidator(round_list(3))] = []
    boltzmann_weights: Annotated[list[float], AfterValidator(round_list(3))] = []

    average_ccs: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    average_ccs_stdev: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
