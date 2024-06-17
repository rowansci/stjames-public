from typing import Any, Optional

from ..mode import Mode
from ..solvent import Solvent
from .workflow import Workflow


class RedoxPotentialWorkflow(Workflow):
    mode: Mode = Mode.RAPID
    solvent: Solvent = Solvent.ACETONITRILE
    reduction: bool = True
    oxidation: bool = True

    # legacy values - remove in future release!
    redox_type: Optional[str] = None
    redox_potential: Optional[float] = None

    # uuids
    neutral_molecule: Optional[str] = None
    anion_molecule: Optional[str] = None
    cation_molecule: Optional[str] = None

    reduction_potential: Optional[float] = None
    oxidation_potential: Optional[float] = None

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
