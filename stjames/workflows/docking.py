"""Docking workflow."""

from typing import Self

from pydantic import BaseModel, ConfigDict, field_validator, model_validator

from ..molecule import Molecule
from ..pdb import PDB
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
    :param initial_molecule: Molecule of interest (currently unused)
    :param mode: Mode for workflow (currently unused)

    New:
    :param molecules: Molecules to dock (optional)
    :param smiles: SMILES strings of the ligands (optional)
    :param target: PDB of the protein
    :param pocket: center (x, y, z) and size (x, y, z) of the pocket

    Results:
    :param scores: docked poses sorted by score
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    molecules: list[Molecule] = []
    smiles: list[str] = []

    target: PDB
    pocket: tuple[Vector3D, Vector3D]

    scores: list[Score] = []

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        """Return a string representation of the Docking workflow."""
        desc = self.target.description
        target = desc.code or desc.title
        ligand = "".join(atom.atomic_symbol for atom in self.initial_molecule.atoms)

        return f"<{type(self).__name__} {target} {ligand}>"

    @model_validator(mode="after")
    def check_molecules(self) -> Self:
        """Check if molecules are provided."""
        if self.molecules and self.smiles:
            raise ValueError("Must provide only one of molecules or smiles, not both")
        elif not self.molecules and not self.smiles:
            raise ValueError("Must provide either molecules or smiles")

        return self

    @field_validator("pocket", mode="after")
    def validate_pocket(cls, pocket: tuple[Vector3D, Vector3D]) -> tuple[Vector3D, Vector3D]:
        center, size = pocket
        if any(q <= 0 for q in size):
            raise ValueError(f"Pocket size must be positive, got: {size}")

        return pocket
