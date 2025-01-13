"""Docking workflow."""

from tempfile import NamedTemporaryFile

import atomium  # type: ignore [import-untyped]
from pydantic import BaseModel, ConfigDict, field_validator

from ..molecule import Molecule
from ..types import Vector3D
from .workflow import Workflow


class Score(BaseModel):
    """
    Pose with its score.

    :param pose: conformation of the ligand when docked
    :param score: score of the pose
    """

    pose: Molecule
    score: float


class DockingWorkflow(Workflow):
    """
    Docking workflow.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param target: pdb of the protein
    :param pocket: center (x, y, z) and size (x, y, z) of the pocket

    Results:
    :param scores: docked poses sorted by score
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    target: atomium.data.File
    pocket: tuple[Vector3D, Vector3D]

    scores: list[Score] = []

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        """Return a string representation of the Docking workflow."""
        target = self.target.title or self.target.name or self.target.id
        ligand = "".join(atom.atomic_symbol for atom in self.initial_molecule.atoms)

        return f"<{type(self).__name__} {target} {ligand}>"

    @field_validator("target", mode="before")
    def read_pdb(cls, v: atomium.data.File | str) -> atomium.data.File:
        """Read the target protein from a PDB file."""
        if isinstance(v, str):
            with NamedTemporaryFile("w", suffix=".pdb") as f:
                f.write(v)
                f.seek(0)
                return atomium.open(f.name)
        elif isinstance(v, atomium.data.File):
            return v

        raise ValueError("Target must be a PDB file or string of pdb")

    @field_validator("pocket", mode="after")
    def validate_pocket(cls, pocket: tuple[Vector3D, Vector3D]) -> tuple[Vector3D, Vector3D]:
        center, size = pocket
        if any(q <= 0 for q in size):
            raise ValueError(f"Pocket size must be positive, got: {size}")

        return pocket
