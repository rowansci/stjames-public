from typing import Annotated

from pydantic import AfterValidator, NonNegativeInt

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
    :param residue_numbers: the indices of the residues on the pocket
    """

    sphere_centers: list[Vector3D]
    sphere_radii: list[float]

    volume: Annotated[float, AfterValidator(round_float(3))]
    score: Annotated[float, AfterValidator(round_float(3))]

    pocket_center: Vector3D
    pocket_sides: Vector3D

    residue_numbers: list[NonNegativeInt]


class PocketDetectionWorkflow(ProteinStructureWorkflow):
    """
    Uses Pocketeer to detect potential binding sites on a protein.

    Inherited:
    :param protein: the protein

    New:
    :param merge_distance: distance for merging pocket spheres, in Å

    Results:
    :param pockets: the located pockets
    """

    merge_distance: Annotated[float, AfterValidator(round_float(3))] = 1.75
    pockets: list[Pocket] = []
