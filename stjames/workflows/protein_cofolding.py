"""Protein cofolding workflow."""

from typing import Annotated, Any, Literal, TypeAlias

from pydantic import AfterValidator, BaseModel, ConfigDict, model_validator

from ..base import LowercaseStrEnum, round_float
from ..types import UUID, round_list
from .workflow import FASTAWorkflow

ProteinUUID: TypeAlias = UUID
CalculationUUID: TypeAlias = UUID


_LDDT_SAMPLE_ROUNDER = round_list(3)


def _round_lddt_samples(values: list[list[float]] | None) -> list[list[float]] | None:
    """Round each diffusion sample's lDDT values to three decimal places."""

    if values is None:
        return None

    return [_LDDT_SAMPLE_ROUNDER(sample) for sample in values]


class CofoldingModel(LowercaseStrEnum):
    """Cofolding model to be used for prediction."""

    CHAI_1R = "chai_1r"
    BOLTZ_1 = "boltz_1"
    BOLTZ_2 = "boltz_2"


class Token(BaseModel):
    """Either an atom in a ligand or a residue in a protein or nucleotide chain."""

    input_type: Literal["ligand", "protein", "dna", "rna"]
    input_index: int
    token_index: int
    atom_name: str | None = None


class ContactConstraint(BaseModel):
    """Contact constraint to be used for prediction."""

    token_1: Token
    token_2: Token
    max_distance: float  # Angstroms
    force: bool = False  # Whether to use potentials to enforce the constraint


class PocketConstraint(BaseModel):
    """Pocket constraint to be used for prediction."""

    input_type: Literal["ligand", "protein", "dna", "rna"]
    input_index: int
    contacts: list[Token]
    max_distance: float  # Angstroms
    force: bool = False  # Whether to use potentials to enforce the constraint


class CofoldingScores(BaseModel):
    """The output scores from co-folding scores."""

    confidence_score: Annotated[float, AfterValidator(round_float(3))]
    ptm: Annotated[float, AfterValidator(round_float(3))]  # predicted template modeling score
    iptm: Annotated[float, AfterValidator(round_float(3))]  # interface predicted template modeling score
    avg_lddt: Annotated[float, AfterValidator(round_float(3))]


class AffinityScore(BaseModel):
    pred_value: Annotated[float, AfterValidator(round_float(3))]
    probability_binary: Annotated[float, AfterValidator(round_float(3))]
    pred_value1: Annotated[float, AfterValidator(round_float(3))]
    probability_binary1: Annotated[float, AfterValidator(round_float(3))]
    pred_value2: Annotated[float, AfterValidator(round_float(3))]
    probability_binary2: Annotated[float, AfterValidator(round_float(3))]


class ProteinCofoldingWorkflow(FASTAWorkflow):
    """
    Workflow for predicting structures.

    Especially protein structures. At least one biological sequence is required.

    Inherited:
    :param initial_protein_sequences: protein sequences of interest
    :param initial_dna_sequences: DNA sequences of interest
    :param initial_rna_sequences: RNA sequences of interest
    :param initial_smiles_list: SMILES strings of interest

    New:
    :param use_msa_server: whether to use the MSA server
    :param use_templates_server: whether to use the templates server
    :param use_potentials: whether to use the potentials (inference-time steering) with Boltz
    :param contact_constraints: Boltz contact constraints
    :param pocket_constraints: Boltz pocket constraints
    :param do_pose_refinement: whether to optimize non-rotatable bonds in output poses
    :param compute_strain: whether to compute the strain of the pose (if `pose_refinement` is enabled)
    :param diffusion_samples: number of samples generated for prediction
    :param model: which cofolding model to use
    :param affinity_score: the affinity score
    :param lddt: per diffusion sample local distance difference test results
    :param predicted_structure_uuid: per diffusion sample UUID of the predicted structure
    :param scores: per diffusion sample cofolding scores
    :param pose: per diffusion sample UUID of the calculation pose
    :param posebusters_valid: per diffusion sample PoseBusters validity flags
    :param strain: per diffusion sample ligand strain, in kcal/mol
    :param predicted_refined_structure_uuid: per diffusion sample UUID of the refined predicted structure (if refinement ran)
    """

    model_config = ConfigDict(validate_assignment=True)

    use_msa_server: bool = False
    use_templates_server: bool = False
    use_potentials: bool = False
    contact_constraints: list[ContactConstraint] = []
    pocket_constraints: list[PocketConstraint] = []
    do_pose_refinement: bool = False
    compute_strain: bool = False
    diffusion_samples: int | None = 1

    model: CofoldingModel = CofoldingModel.BOLTZ_2
    affinity_score: AffinityScore | None = None
    lddt: Annotated[list[list[float]] | None, AfterValidator(_round_lddt_samples)] = None

    predicted_structure_uuid: list[ProteinUUID] | None = None
    scores: list[CofoldingScores] | None = None
    pose: list[CalculationUUID] | None = None
    posebusters_valid: list[bool] | None = None
    strain: Annotated[list[float] | None, AfterValidator(round_list(3))] = None
    predicted_refined_structure_uuid: list[ProteinUUID] | None = None

    @model_validator(mode="before")
    @classmethod
    def _ensure_list_outputs(cls, values: Any) -> Any:
        """Allow single output values by wrapping them in a list."""

        if not isinstance(values, dict):
            return values

        list_fields = (
            "predicted_structure_uuid",
            "scores",
            "pose",
            "posebusters_valid",
            "strain",
            "predicted_refined_structure_uuid",
        )

        for field in list_fields:
            value = values.get(field)
            if value is None or isinstance(value, list):
                continue
            values[field] = [value]

        lddt_values = values.get("lddt")
        if lddt_values is None:
            return values

        if not isinstance(lddt_values, list):
            values["lddt"] = [lddt_values]
            return values

        if not lddt_values:
            return values

        if all(isinstance(sample, list) for sample in lddt_values):
            return values

        values["lddt"] = [lddt_values]
        return values
