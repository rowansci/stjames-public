from typing import Annotated

from pydantic import AfterValidator

from ..base import round_float, round_optional_float
from ..conformers import ConformerGenSettingsUnion, ETKDGSettings
from ..mode import Mode
from ..settings import Settings
from ..types import UUID
from .multistage_opt import MultiStageOptSettings
from .workflow import MoleculeWorkflow


class StrainWorkflow(MoleculeWorkflow):
    """
    Workflow for calculating the strain of a given molecular geometry.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param conf_gen_settings: conformer-search settings
    :param multistage_opt_settings: optimization settings
    :param harmonic_constraint_spring_constant: spring constant for constraints, in kcal/mol/Å
    :param constrain_hydrogens: whether or not to constrain hydrogens

    Results:
    :param conformers: list of conformer UUIDs
    :param constrained_optimization: UUID of optimized strained structure
    :param strain: actual strain, in kcal/mol
    """

    conf_gen_settings: ConformerGenSettingsUnion = ETKDGSettings(max_confs=50)
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
