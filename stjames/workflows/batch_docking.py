"""High-throughput docking workflow."""

from typing import Annotated, Any

from pydantic import AfterValidator, ConfigDict, field_validator, model_validator

from ..pdb import PDB
from ..types import UUID, Vector3D, round_list
from .docking import VinaSettings
from .workflow import BatchSMILESWorkflow, ProteinStructureWorkflow


class BatchDockingWorkflow(BatchSMILESWorkflow, ProteinStructureWorkflow):
    """
    Batch docking workflow.

    Inherited:
    :param initial_smiles_list: list of SMILES
    :param protein: PDB of the protein, or the UUID of the protein.

    New:
    :param target: PDB of the protein, or the UUID of the protein. DEPRECATED
    :param pocket: center (x, y, z) and size (x, y, z) of the pocket
    :param docking_settings: how to run each docking calculation

    Results:
    :param best_scores: best score for each SMILES string
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    target: PDB | UUID
    pocket: tuple[Vector3D, Vector3D]

    docking_settings: VinaSettings = VinaSettings()
    best_scores: Annotated[list[float | None], AfterValidator(round_list(3))] = []

    @model_validator(mode="before")
    def harmonize_target_and_protein(cls, data: Any) -> Any:  # noqa: N805
        """
        Syncs data between "target" and "protein" field.
        """
        protein = data.get("protein")
        target = data.get("target")

        if target and not protein:
            data["protein"] = target
        elif not target and protein:
            data["target"] = protein

        return data

    @field_validator("pocket", mode="after")
    def validate_pocket(cls, pocket: tuple[Vector3D, Vector3D]) -> tuple[Vector3D, Vector3D]:
        _center, size = pocket
        if any(q <= 0 for q in size):
            raise ValueError(f"Pocket size must be positive, got: {size}")
        return pocket
