from typing import Annotated

from pydantic import AfterValidator

from ..base import Base, round_float
from ..conformers import ConformerClusteringSettings, ConformerGenSettingsUnion, ConformerProperties, ConformerSearchMixin, iMTDSettings
from ..method import Method
from ..mode import Mode
from ..settings import Settings
from ..solvent import Solvent, SolventModel, SolventSettings
from ..types import UUID
from .multistage_opt import MultiStageOptSettings
from .workflow import MoleculeWorkflow


class SolventDependentConformer(Base):
    """
    Stores a single conformer scored in many different solvents.

    :param calculation: conformer (as a calculation)
    :param free_energy_by_solvent: free energy in every solvent, in Hartree
    :param relative_free_energy_by_solvent: relative free energy vs. lowest-energy conformer in every solvent, in kcal/mol
    :param population_by_solvent: population in every solvent (number between 0 and 1)
    """

    calculation: UUID
    free_energy_by_solvent: dict[Solvent, Annotated[float, AfterValidator(round_float(3))]]
    relative_free_energy_by_solvent: dict[Solvent, Annotated[float, AfterValidator(round_float(3))]]
    population_by_solvent: dict[Solvent, Annotated[float, AfterValidator(round_float(3))]]


class SolventDependentConformersWorkflow(ConformerSearchMixin, MoleculeWorkflow):
    """
    Conformers are generated through `conf_gen_settings`.

    Clustering is then performed, and optimization is conducted.
    The optimization of the conformers is governed by `multistage_opt_settings`.

    Final conformer scoring is done through a multi-level scheme:
        - A single-point gas-phase energy run through `multistage_opt_settings.sp_settings`.
        - A separate thermal free-energy correction computed via a single-point Hessian (GFN2-xTB).
        - A per-solvent CPCM-X calculation.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param conf_gen_settings: settings for conformer generation
    :param multistage_opt_settings: set by mode unless mode=MANUAL (ignores additional settings if set)

    New:
    :param solvents: solvents to study
    :param conformers: output conformers with per-solvent energies and weights
    :param per_solvent_properties: metrics for how overall distribution changes by solvent
    :param relative_free_energy_by_solvent: how free energy changes by solvent, for predicting ∆G_transfer
    """

    solvents: list[Solvent] = [
        Solvent.HEXANE,
        Solvent.OCTANOL,
        Solvent.CHLOROFORM,
        Solvent.DIMETHYLSULFOXIDE,
        Solvent.WATER,
    ]

    conf_gen_settings: ConformerGenSettingsUnion | None = iMTDSettings(
        max_confs=None,
        speed="normal",
        solvent_settings=SolventSettings(solvent=Solvent.WATER, model=SolventModel.ALPB),
        reopt=False,
        mtd_method=Method.GFN2_XTB,
    )

    conformer_clustering_settings: ConformerClusteringSettings | None = ConformerClusteringSettings()

    multistage_opt_settings: MultiStageOptSettings = MultiStageOptSettings(
        optimization_settings=[
            Settings(
                method="gfn2_xtb",
                solvent_settings=SolventSettings(
                    solvent=Solvent.WATER,
                    model=SolventModel.ALPB,
                ),
                tasks=["optimize"],
            )
        ],
        sp_settings=Settings(method="g_xtb", tasks=["energy"]),
        mode=Mode.MANUAL,
    )

    conformers: list[SolventDependentConformer] = []
    per_solvent_properties: dict[Solvent, ConformerProperties] = {}
    relative_free_energy_by_solvent: dict[Solvent, Annotated[float, AfterValidator(round_float(3))]] = {}
