"""Protein Cofolding Workflow."""

from typing import Annotated, Literal

from pydantic import AfterValidator, BaseModel, ConfigDict

from ..base import LowercaseStrEnum, round_float
from ..types import UUID, round_optional_list
from .workflow import FASTAWorkflow


class CofoldingModel(LowercaseStrEnum):
    """Cofolding model to be used for prediction."""

    CHAI_1R = "chai_1r"
    BOLTZ_1 = "boltz_1"
    BOLTZ_2 = "boltz_2"


class Token(BaseModel):
    """Either a atom in a ligand or a residue in a protein."""

    input_type: Literal["ligand", "protein"]
    input_index: int
    token_index: int


class ContactConstraint(BaseModel):
    """Contact constraint to be used for prediction."""

    token_1: Token
    token_2: Token
    max_distance: float  # Angstroms
    force: bool = False  # Whether to use potentials to enforce the constraint


class PocketConstraint(BaseModel):
    """Pocket constraint to be used for prediction."""

    input_type: Literal["ligand", "protein"]
    input_index: int
    contacts: list[Token]
    max_distance: float  # Angstroms
    force: bool = False  # Whether to use potentials to enforce the constraint


class CofoldingScores(BaseModel):
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
    A workflow for predicting structures. Especially protein structures.

    Inherited:
    :param initial_protein_sequences: protein sequences of interest
    :param initial_smiles_list: SMILES strings of interest

    New:
    :param use_msa_server: whether to use the MSA server
    :param use_templates_server: whether to use the templates server
    :param predicted_structure_uuid: UUID of the predicted structure
    """

    model_config = ConfigDict(validate_assignment=True)

    use_msa_server: bool = False
    use_templates_server: bool = False
    use_potentials: bool = False
    contact_constraints: list[ContactConstraint] = []
    pocket_constraints: list[PocketConstraint] = []

    predicted_structure_uuid: UUID | None = None
    scores: CofoldingScores | None = None
    model: CofoldingModel = CofoldingModel.BOLTZ_2
    affinity_score: AffinityScore | None = None
    lddt: Annotated[list[float] | None, AfterValidator(round_optional_list(3))] = None
