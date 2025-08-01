"""Nuclear-magnetic-resonance spectroscopy workflow."""

from typing import Annotated

from pydantic import AfterValidator

from ..base import LowercaseStrEnum
from ..mode import Mode
from ..settings import Settings
from ..solvent import Solvent
from ..types import UUID, round_list
from .conformer_search import ConformerGenSettings, iMTDSettings
from .multistage_opt import MultiStageOptSettings
from .workflow import MoleculeWorkflow


class NMRMethod(LowercaseStrEnum):
    MAGNETZERO = "magnet-zero"


class NMRSpectroscopyWorkflow(MoleculeWorkflow):
    """
    Workflow for calculating NMR spectra.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param nmr_method: how to run the NMR calculations
    :param solvent: the solvent in which to run the calculations
    :param conf_gen_settings : the conformer-search settings. if `None`, no conformer search will be performed
    :param multistage_opt_settings: the optimization settings. if `None`, no optimization will be performed

    Results:
    :param conformers: list of conformer UUIDs
    :param boltzmann_weights: the boltzmann weights for each conformer
    :param per_conformer_isotopic_shieldings: the per-atom shieldings for each conformer
    :param isotopic_shieldings: the per-atom shieldings
    """

    nmr_method: NMRMethod = NMRMethod.MAGNETZERO
    solvent: Solvent = Solvent.CHLOROFORM

    conf_gen_settings: ConformerGenSettings | None = iMTDSettings(mode="careful")
    multistage_opt_settings: MultiStageOptSettings | None = MultiStageOptSettings(
        mode=Mode.MANUAL,
        optimization_settings=[Settings(method="aimnet2_wb97md3")],
    )

    conformers: list[UUID] = []
    boltzmann_weights: Annotated[list[float], AfterValidator(round_list(3))] = []
    per_conformer_isotropic_shieldings: list[Annotated[list[float], AfterValidator(round_list(3))]] = []
    isotropic_shieldings: Annotated[list[float], AfterValidator(round_list(3))] = []
