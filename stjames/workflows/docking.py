"""Docking workflow."""

from typing import Annotated, Self, TypeAlias

from pydantic import AfterValidator, ConfigDict, field_validator, model_validator

from ..base import Base, round_float
from ..pdb import PDB
from ..types import UUID, Vector3D
from .conformer_search import ConformerGenSettingsUnion, ETKDGSettings
from .workflow import MoleculeWorkflow

ProteinUUID: TypeAlias = UUID
CalculationUUID: TypeAlias = UUID


class Score(Base):
    """
    Pose with its score.

    :param pose: conformation of the ligand when docked (calculation UUID)
    :param complex_pdb: the UUID of the protein–ligand complex (protein UUID)
    :param score: score of the pose, in kcal/mol
    :param posebusters_valid: whether or not the ligand pose passes the PoseBusters tests
    :param strain: strain in kcal/mol
    """

    pose: CalculationUUID | None
    complex_pdb: ProteinUUID | None
    score: Annotated[float, AfterValidator(round_float(3))]
    posebusters_valid: bool
    strain: float | None


class DockingSettings(Base):
    """
    Base class for controlling how docked poses are generated.

    :param max_poses: the maximum number of poses generated per input molecule
    """

    max_poses: int = 4


class VinaSettings(DockingSettings):
    """
    Controls how AutoDock Vina is run.

    :param exhaustiveness: how many times Vina attempts to find a pose.
        8 is typical, 32 is considered relatively careful.
    """

    exhaustiveness: int = 8


class DockingWorkflow(MoleculeWorkflow):
    """
    Docking workflow.

    Note that the protein can be supplied either by UUID or raw PDB object.
    We anticipate that the former will dominate deployed usage, but the latter is handy for isolated testing.
    If, for whatever reason, the workflow is initialized with both a `target_uuid` and a `target`, the UUID will be ignored.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param do_csearch: whether to csearch starting structures
    :param conformer_gen_settings: settings for initial conformer search.
    :param do_optimization: whether to optimize starting structures
    :param optimization_settings: settings for conformer optimization.
    :param do_pose_refinement: whether to optimize non-rotatable bonds in output poses
    :param target: PDB of the protein.
    :param target_uuid: UUID of the protein.
    :param pocket: center (x, y, z) and size (x, y, z) of the pocket

    Results:
    :param conformers: UUIDs of optimized conformers
    :param scores: docked poses sorted by score
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    target: PDB | None = None
    target_uuid: UUID | None = None
    pocket: tuple[Vector3D, Vector3D]

    docking_settings: VinaSettings = VinaSettings()

    do_csearch: bool = True
    conformer_gen_settings: ConformerGenSettingsUnion = ETKDGSettings(mode="reckless")
    do_optimization: bool = True
    do_pose_refinement: bool = True

    conformers: list[CalculationUUID] = []
    scores: list[Score] = []

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        """Return a string representation of the Docking workflow."""
        if self.target is not None:
            desc = self.target.description
            target = desc.code or desc.title
        else:
            target = ""

        ligand = "".join(atom.atomic_symbol for atom in self.initial_molecule.atoms)
        return f"<{type(self).__name__} {target} {ligand}>"

    @model_validator(mode="after")
    def check_protein(self) -> Self:
        """Check if protein is provided."""
        if not self.target and not self.target_uuid:
            raise ValueError("Must provide either target or target_uuid")
        return self

    @field_validator("pocket", mode="after")
    def validate_pocket(cls, pocket: tuple[Vector3D, Vector3D]) -> tuple[Vector3D, Vector3D]:
        _center, size = pocket
        if any(q <= 0 for q in size):
            raise ValueError(f"Pocket size must be positive, got: {size}")
        return pocket
