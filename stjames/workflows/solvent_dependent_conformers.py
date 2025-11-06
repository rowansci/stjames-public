from pydantic import PositiveFloat

from ..mode import Mode
from ..settings import Settings
from ..solvent import Solvent, SolventModel, SolventSettings
from ..types import UUID
from .conformer_search import ConformerClusteringSettings, ConformerSearchMixin, MonteCarloMultipleMinimumSettings
from .multistage_opt import MultiStageOptMixin, MultiStageOptSettings
from .workflow import MoleculeWorkflow


class SolventDependentConformer:
    molecule: UUID
    free_energy_by_solvent: dict[Solvent, float]
    relative_free_energy_by_solvent: dict[Solvent, float]
    population_by_solvent: dict[Solvent, PositiveFloat]


class ConformerEnsembleProperties:
    mean_solvent_accessible_surface_area: PositiveFloat
    mean_polar_solvent_accessible_surface_area: PositiveFloat
    mean_radius_of_gyration: PositiveFloat


class SolventDependentConformersWorkflow(MultiStageOptMixin, ConformerSearchMixin, MoleculeWorkflow):
    """

    The optimization of the conformers is governed by `multistage_opt_settings`.

    Final conformer scoring is done through a multi-level scheme:
        - A single-point gas-phase energy run through `multistage_opt_settings.sp_settings`.
        - A separate thermal free-energy correcting computed via a single-point Hessian (GFN2-xTB).
        - A per-solvent CPCM-X calculation.



    (more here)

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

    conf_gen_mode = Mode.MANUAL
    conf_gen_settings = MonteCarloMultipleMinimumSettings(num_iterations=500)

    conformer_clustering_settings = ConformerClusteringSettings()

    multistage_opt_settings: MultiStageOptSettings = MultiStageOptSettings(
        optimization_settings=[
            Settings(
                method="gfn2_xtb",
                solvent_settings=SolventSettings(
                    solvent=Solvent.WATER,
                    model=SolventModel.ALPB,
                ),
            )
        ],
        sp_settings=Settings(method="g_xtb"),
    )

    conformers: list[SolventDependentConformer]
    per_solvent_properties: dict[Solvent, ConformerEnsembleProperties]
    relative_free_energy_by_solvent: dict[Solvent, float]
