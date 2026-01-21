"""Docking workflow."""

from typing import Annotated, Any, Literal, Self, TypeAlias

from pydantic import AfterValidator, ConfigDict, field_validator, model_validator

from ..base import Base, round_float
from ..conformers import ConformerGenSettingsUnion, ETKDGSettings
from ..pdb import PDB
from ..types import UUID, Vector3D
from .workflow import MoleculeWorkflow, ProteinStructureWorkflow

ProteinUUID: TypeAlias = UUID
CalculationUUID: TypeAlias = UUID


class Score(Base):
    """
    Pose with its score.

    :param pose: conformation of the ligand when docked (calculation UUID)
    :param complex_pdb: UUID of protein–ligand complex (protein UUID)
    :param score: score of the pose, (kcal/mol)
    :param posebusters_valid: whether or not the ligand pose passes the PoseBusters tests
    :param strain: strain (kcal/mol)
    :param rmsd: RMSD from the reference, if there's a reference molecule to dock against (Å)
    """

    pose: CalculationUUID | None
    complex_pdb: ProteinUUID | None
    score: Annotated[float, AfterValidator(round_float(3))]
    posebusters_valid: bool
    strain: float | None
    rmsd: Annotated[float, AfterValidator(round_float(3))] | None = None


class DockingSettings(Base):
    """
    Base class for controlling how docked poses are generated.

    :param max_poses: maximum number of poses generated per input molecule
    """

    max_poses: int = 4


class VinaSettings(DockingSettings):
    """
    Controls how AutoDock Vina is run.

    :param executable: which Vina implementation is run.
    :param scoring_function: which scoring function is employed.
    :param exhaustiveness: how many times Vina attempts to find a pose.
        8 is typical, 32 is considered relatively careful.
    """

    executable: Literal["qvina2", "qvina-w", "vina"] = "vina"
    scoring_function: Literal["vinardo", "vina"] = "vinardo"
    exhaustiveness: int = 8

    @model_validator(mode="after")
    def check_executable_scoring_function(self) -> Self:
        """Check if the combination of exectuable and scoring function is supported."""
        if (self.executable in {"qvina2", "qvina-w"}) and (self.scoring_function == "vinardo"):
            raise ValueError("QVina does not implement the Vinardo scoring function!")
        return self


class DockingWorkflow(MoleculeWorkflow, ProteinStructureWorkflow):
    """
    Docking workflow.

    Note that the protein can be supplied either by UUID or raw PDB object.
    We anticipate that the former will dominate deployed usage, but the latter is handy for isolated testing.
    If, for whatever reason, the workflow is initialized with both a `target_uuid` and a `target`, the UUID will be ignored.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)
    :param protein: Protein target, as PDB or UUID

    New:
    :param target: PDB of the protein; DEPRECATED.
    :param target_uuid: UUID of the protein; DEPRECATED.
    :param pocket: center (x, y, z) and size (x, y, z) of the pocket
    :param do_csearch: whether to csearch starting structures
    :param conformer_gen_settings: settings for initial conformer search.
    :param do_optimization: whether to optimize starting structures
    :param do_pose_refinement: whether to optimize non-rotatable bonds in output poses

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
    conformer_gen_settings: ConformerGenSettingsUnion = ETKDGSettings(
        num_initial_confs=200,
        num_confs_considered=50,
        max_confs=20,
        max_mmff_energy=20,
    )
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

    @model_validator(mode="before")
    def harmonize_target_and_protein(cls, data: Any) -> Any:  # noqa: N805
        """
        Syncs data between "target"/"target_uuid" and "protein" field.
        """
        protein = data.get("protein", False)
        target = data.get("target", False)
        target_uuid = data.get("target_uuid", False)

        if not protein:
            if target:
                data["protein"] = target
            elif target_uuid:
                data["protein"] = target_uuid
        elif not target and not target_uuid:
            if isinstance(data["protein"], PDB):
                data["target"] = protein
            elif isinstance(data["protein"], UUID):
                data["target_uuid"] = protein

        return data

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


class AnalogueDockingWorkflow(MoleculeWorkflow, ProteinStructureWorkflow):
    """
    Workflow for docking analogues:
    (1) Conformers are generated in analogous poses to the initial molecule.
    (2) They're then optimized locally using the docking scoring function.
    (3) PoseBusters is used to check the validity of the output poses.

    Inherited:
    :param initial_molecule: molecule of interest, to which subsequent molecules will be aligned
    :param mode: Mode for workflow (currently unused)
    :param protein: PDB or UUID

    New:
    :param analogues: SMILES for analogues of `initial_molecule`
    :param docking_settings: how docking should be run

    Results:
    :param analogue_scores: docked poses for each analogue of form {smiles: list[poses]}
    """

    analogues: list[str]
    docking_settings: VinaSettings = VinaSettings()

    analogue_scores: dict[str, list[Score]] = {}
