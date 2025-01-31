"""Docking workflow."""

from typing import Annotated

from pydantic import AfterValidator, ConfigDict, field_validator

from ..base import Base, round_float
from ..pdb import PDB
from ..types import UUID, Vector3D
from .workflow import Workflow


class Score(Base):
    """
    Pose with its score.

    :param pose: conformation of the ligand when docked
    :param score: score of the pose, in kcal/mol
    """

    pose: UUID | None  # for calculation
    score: Annotated[float, AfterValidator(round_float(3))]


class DockingWorkflow(Workflow):
    """
    Docking workflow.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param molecules: Molecules to dock (optional)
    :param smiles: SMILES strings of the ligands (optional)
    :param do_csearch: whether to csearch starting structures
    :param do_optimization: whether to optimize starting structures
    :param conformers: UUIDs of optimized conformers
    :param target: PDB of the protein
    :param pocket: center (x, y, z) and size (x, y, z) of the pocket

    Results:
    :param scores: docked poses sorted by score
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    do_csearch: bool = True
    do_optimization: bool = True
    conformers: list[UUID] = []

    target: PDB
    pocket: tuple[Vector3D, Vector3D]

    do_pose_hydrogen_refinement: bool = True
    scores: list[Score] = []

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        """Return a string representation of the Docking workflow."""
        desc = self.target.description
        target = desc.code or desc.title
        ligand = "".join(atom.atomic_symbol for atom in self.initial_molecule.atoms)

        return f"<{type(self).__name__} {target} {ligand}>"

    @field_validator("pocket", mode="after")
    def validate_pocket(cls, pocket: tuple[Vector3D, Vector3D]) -> tuple[Vector3D, Vector3D]:
        center, size = pocket
        if any(q <= 0 for q in size):
            raise ValueError(f"Pocket size must be positive, got: {size}")

        return pocket
