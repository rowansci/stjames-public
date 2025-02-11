"""ADME-Tox property prediction workflow."""

from .workflow import MoleculeWorkflow


class ADMETWorkflow(MoleculeWorkflow):
    """
    A workflow for predicting ADME-Tox properties.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param properties: predicted properties
    """

    properties: dict[str, float | int] | None = None
