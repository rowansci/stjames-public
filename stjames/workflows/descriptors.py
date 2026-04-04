"""Molecular descriptors workflow."""

from ..solvent import Solvent
from ..types import UUID
from .workflow import MoleculeWorkflow

Descriptors = dict[str, dict[str, float] | tuple[float | None, ...] | float]


class DescriptorsWorkflow(MoleculeWorkflow):
    """
    A workflow for calculating molecular descriptors.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (not used)

    New:
    :param do_optimization: whether to optimize with GFN2-xTB before calculating descriptors
    :param solvent: solvent to use for optimization and descriptor calculation
    :param optimization: UUID of optimization
    :param descriptors: calculated descriptors
    """

    do_optimization: bool = True
    solvent: Solvent | None = None

    optimization: UUID | None = None
    descriptors: Descriptors | None = None
