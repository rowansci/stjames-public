from typing import Any

from ..types import UUID
from .multistage_opt import MultiStageOptMixin
from .workflow import Workflow


class RedoxPotentialWorkflow(Workflow, MultiStageOptMixin):
    """
    Workflow for computing spin states of molecules.

    Uses the modes from MultiStageOptSettings.

    Influenced by:
    [Performance of Quantum Chemistry Methods for Benchmark Set of Spin–State
    Energetics Derived from Experimental Data of 17 Transition Metal Complexes
    (SSE17)](https://chemrxiv.org/engage/chemrxiv/article-details/66a8b15cc9c6a5c07a792487)

    Inherited
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow
    :param multistage_opt_settings: set by mode unless mode=MANUAL (ignores additional settings if set)
    :param solvent: solvent to use for optimization
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param constraints: constraints to add
    :param transition_state: whether this is a transition state
    :param frequencies: whether to calculate frequencies

    Overridden:
    :param mso_mode: Mode for MultiStageOptSettings

    New:
    :param reduction: whether or not to calculate the reduction half-reaction
    :param oxidation: whether or not to calculate the oxidation half-reaction
    :param neutral_molecule: UUID of the calculation for the neutral molecule
    :param anion_molecule: UUID of the calculation for the anion molecule
    :param cation_molecule: UUID of the calculation for the cation molecule
    :param reduction_potential: the final potential, in V
    :param oxidation_potential: the final potential, in V

    Legacy:
    :param redox_type: one of "reduction" or "oxidation"
    :param redox_potential: the corresponding potential, in V
    """

    reduction: bool = True
    oxidation: bool = True

    # legacy values - remove in future release!
    redox_type: str | None = None
    redox_potential: float | None = None

    # uuids
    neutral_molecule: UUID | None = None
    anion_molecule: UUID | None = None
    cation_molecule: UUID | None = None

    reduction_potential: float | None = None
    oxidation_potential: float | None = None

    def model_post_init(self, __context: Any) -> None:
        """Keep back-compatible with old schema."""
        if self.redox_type == "oxidation":
            self.oxidation = True
            self.reduction = False
            if self.oxidation_potential is None:
                self.oxidation_potential = self.redox_potential
        elif self.redox_type == "reduction":
            self.oxidation = False
            self.reduction = True
            if self.reduction_potential is None:
                self.reduction_potential = self.redox_potential
