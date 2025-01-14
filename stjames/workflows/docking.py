"""Docking workflow."""

from tempfile import NamedTemporaryFile
from typing import Self

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

    target: PDBMolecule
    pocket: tuple[Vector3D, Vector3D]

    scores: list[Score] = []

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        """Return a string representation of the Docking workflow."""
        target = self.target.title or self.target.name or self.target.id
        ligand = "".join(atom.atomic_symbol for atom in self.initial_molecule.atoms)

        return f"<{type(self).__name__} {target} {ligand}>"

    @model_validator
    @classmethod
    def check_molecules(self) -> Self:
        """Check if molecules are provided."""
        if self.molecules and self.smiles:
            raise ValueError("Must provide only one of molecules or smiles, not both")
        elif not self.molecules and not self.smiles:
            raise ValueError("Must provide either molecules or smiles")

        return self

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
