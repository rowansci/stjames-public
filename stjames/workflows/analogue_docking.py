"""Analogue docking workflow."""

from pydantic import ConfigDict

from .docking import Score, VinaSettings
from .workflow import MoleculeWorkflow, ProteinStructureWorkflow


class AnalogueDockingWorkflow(MoleculeWorkflow, ProteinStructureWorkflow):
    """
    Workflow for docking analogues:
    (1) Conformers are generated in analogous poses to the initial molecule.
    (2) They're then optimized locally using the docking scoring function.
    (3) PoseBusters is used to check the validity of the output poses.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param protein: PDB or UUID
    :param mode: Mode for workflow (currently unused)

    New:
    :param analogues: the SMILES for the analogues
    :param docking_settings: how docking should be run

    Results:
    :param analogue_scores: the docked poses for each analogue
        (The SMILES string from above is the key, and the list of poses is the value.)
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    analogues: list[str]
    docking_settings: VinaSettings = VinaSettings()

    analogue_scores: dict[str, list[Score]] = {}
