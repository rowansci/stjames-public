"""Nuclear-magnetic-resonance spectroscopy workflow."""

from typing import Annotated

from pydantic import AfterValidator

from ..base import Base, LowercaseStrEnum, round_float
from ..mode import Mode
from ..settings import Settings
from ..solvent import Solvent
from ..types import UUID, round_list
from .conformer_search import ConformerGenSettingsUnion, iMTDSettings
from .multistage_opt import MultiStageOptSettings
from .workflow import MoleculeWorkflow


class NMRMethod(LowercaseStrEnum):
    MAGNETZERO = "magnet-zero"


class NMRPeak(Base):
    """
    Represents a single NMR peak.

    :param nucleus: the atomic number of the nucleus in question
    :param shift: the chemical shift of the peak
    :param atom_indices: the zero-indices of the atoms giving rise to the peak
    """

    nucleus: int
    shift: Annotated[float, AfterValidator(round_float(3))]
    atom_indices: list[int]


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
    :param per_conformer_chemical_shifts: the per-atom shifts for each conformer
    :param chemical_shifts: the per-atom shifts
    :param symmetry_equivalent_nuclei: 0-indexed atoms which are equivalent to one another
    :param predicted_peaks: the predicted NMR peaks
    """

    nmr_method: NMRMethod = NMRMethod.MAGNETZERO
    solvent: Solvent = Solvent.CHLOROFORM

    conf_gen_settings: ConformerGenSettingsUnion | None = iMTDSettings(mode="careful")
    multistage_opt_settings: MultiStageOptSettings | None = MultiStageOptSettings(
        mode=Mode.MANUAL,
        optimization_settings=[Settings(method="aimnet2_wb97md3", tasks=["optimize"])],
    )

    conformers: list[UUID] = []
    boltzmann_weights: Annotated[list[float], AfterValidator(round_list(3))] = []
    per_conformer_chemical_shifts: list[Annotated[list[float | None], AfterValidator(round_list(3))]] = []
    chemical_shifts: Annotated[list[float | None], AfterValidator(round_list(3))] = []
    symmetry_equivalent_nuclei: list[list[int]] = []

    predicted_peaks: dict[int, list[NMRPeak]] = {}
