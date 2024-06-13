from typing import Optional

from ..solvent import Solvent
from .workflow import Workflow


class RedoxPotentialWorkflow(Workflow):
    solvent: Solvent = Solvent.ACETONITRILE
    reduction: bool = True
    oxidation: bool = True

    # uuids
    neutral_molecule: Optional[str] = None
    anion_molecule: Optional[str] = None
    cation_molecule: Optional[str] = None

    reduction_potential: Optional[float] = None
    oxidation_potential: Optional[float] = None
