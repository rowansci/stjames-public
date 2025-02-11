"""Molecular descriptors workflow."""

from ..types import UUID
from .workflow import MoleculeWorkflow

Descriptors = dict[str, dict[str, float] | tuple[float | None, ...] | float]


class DescriptorsWorkflow(MoleculeWorkflow):
    """
    A workflow for calculating molecular descriptors.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow

    New:
    :param optimization: UUID of optimization
    :param descriptors: calculated descriptors
    """

    optimization: UUID | None = None

    descriptors: Descriptors | None = None
