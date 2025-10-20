from typing import Annotated

from pydantic import AfterValidator

from ..base import round_float, round_optional_float
from ..mode import Mode
from ..settings import Settings
from ..types import UUID
from .conformer_search import ConformerGenSettingsUnion, ETKDGSettings
from .multistage_opt import MultiStageOptSettings
from .workflow import MoleculeWorkflow


class StrainWorkflow(MoleculeWorkflow):
    """
    Workflow for calculating the strain of a given molecular geometry.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param conf_gen_settings : the conformer-search settings.
    :param multistage_opt_settings: the optimization settings.
    :param harmonic_constraint_spring_constant: the spring constant for the constraints in kcal/mol/Å
    :param constrain_hydrogens: whether or not to constrain hydrogens

    Results:
    :param conformers: list of conformer UUIDs
    :param constrained_optimization: the UUID of the optimized strained structure
    :param strain: the actual strain in kcal/mol
    """

    conf_gen_settings: ConformerGenSettingsUnion = ETKDGSettings(mode="rapid")
    multistage_opt_settings: MultiStageOptSettings = MultiStageOptSettings(
        mode=Mode.MANUAL,
        optimization_settings=[Settings(method="aimnet2_wb97md3", tasks=["optimize"])],
        singlepoint_settings=Settings(method="aimnet2_wb97md3", tasks=["energy"], solvent_settings={"solvent": "water", "model": "cpcmx"}),
    )

    harmonic_constraint_spring_constant: Annotated[float, AfterValidator(round_float(3))] = 5.0
    constrain_hydrogens: bool = False

    constrained_optimization: UUID | None = None
    conformers: list[UUID | None] = []
    strain: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
