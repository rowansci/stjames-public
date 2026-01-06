from typing import Annotated, Self

from pydantic import AfterValidator, model_validator

from ..base import Base, round_float
from ..types import Vector3D
from .workflow import ProteinStructureWorkflow


class Pocket(Base):
    """
    Represents a pocket.

    :param sphere_centers: the centers of the detected spheres
    :param sphere_radii: the radii of the detected spheres
    :param volume: the volume, in Å**3
    :param score: the druggability / quality score, larger scores are better
    :param pocket_center: the center of the bounding box
    :param pocket_sides: the side lengths of the bounding box
    """

    sphere_centers: list[Vector3D]
    sphere_radii: list[float]

    volume: Annotated[float, AfterValidator(round_float(3))] = 1.75
    score: Annotated[float, AfterValidator(round_float(3))] = 1.75

    pocket_center: list[Vector3D]
    pocket_sides: list[Vector3D]


class PocketDetectionWorkflow(ProteinStructureWorkflow):
    """
    Uses Pocketeer to detect potential binding sites on a protein.

    Inherited:
    :param protein: the protein
    :param protein_uuid: the protein's UUID

    New:
    :param merge_distance: distance for merging pocket spheres, in Å

    Results:
    :param pockets: the located pockets
    """

    merge_distance: Annotated[float, AfterValidator(round_float(3))] = 1.75
    pockets: list[Pocket] = []

    @model_validator(mode="after")
    def check_protein(self) -> Self:
        """Check if protein is provided."""
        if not self.protein and not self.protein_uuid:
            raise ValueError("Must provide either target or target_uuid")
        return self
