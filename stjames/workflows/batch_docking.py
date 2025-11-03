"""High-throughput docking workflow."""

from typing import Annotated

from pydantic import AfterValidator, ConfigDict, field_validator

from ..pdb import PDB
from ..types import UUID, Vector3D, round_list
from .docking import VinaSettings
from .workflow import BatchSMILESWorkflow


class BatchDockingWorkflow(BatchSMILESWorkflow):
    """
    Docking workflow.

    Note that the protein can be supplied either by UUID or raw PDB object.
    We anticipate that the former will dominate deployed usage, but the latter is handy for isolated testing.
    If, for whatever reason, the workflow is initialized with both a `target_uuid` and a `target`, the UUID will be ignored.

    Inherited:
    :param initial_smiles_list: list of SMILES

    New:
    :param target: PDB of the protein, or the UUID of the protein.
    :param pocket: center (x, y, z) and size (x, y, z) of the pocket
    :param docking_settings: how to run each docking calculation

    Results:
    :param best_scores: the best score for each SMILES string
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    target: PDB | UUID
    pocket: tuple[Vector3D, Vector3D]

    docking_settings: VinaSettings = VinaSettings()
    best_scores: Annotated[list[float | None], AfterValidator(round_list(3))] = []

    @field_validator("pocket", mode="after")
    def validate_pocket(cls, pocket: tuple[Vector3D, Vector3D]) -> tuple[Vector3D, Vector3D]:
        _center, size = pocket
        if any(q <= 0 for q in size):
            raise ValueError(f"Pocket size must be positive, got: {size}")
        return pocket
