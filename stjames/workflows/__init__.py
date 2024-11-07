# ruff: noqa: F405
from typing import Literal

from .admet import *
from .basic_calculation import *
from .bde import *
from .conformer import *
from .conformer_search import *
from .descriptors import *
from .fukui import *
from .molecular_dynamics import *
from .multistage_opt import *
from .pka import *
from .redox_potential import *
from .scan import *
from .spin_states import *
from .tautomer import *

WORKFLOW_NAME = Literal[
    "admet",
    "basic_calculation",
    "bde",
    "conformers",
    "conformer_search",
    "descriptors",
    "fukui_index",
    "molecular_dynamics",
    "multistage_opt",
    "pka",
    "redox_potential",
    "scan",
    "spin_states",
    "tautomers",
]

WORKFLOW_MAPPING = {
    "admet": ADMETWorkflow,
    "basic_calculation": BasicCalculationWorkflow,
    "bde": BDEWorkflow,
    "conformers": ConformerWorkflow,
    "conformer_search": ConformerSearchWorkflow,
    "descriptors": DescriptorsWorkflow,
    "fukui_index": FukuiIndexWorkflow,
    "molecular_dynamics": MolecularDynamicsWorkflow,
    "multistage_opt": MultiStageOptWorkflow,
    "pka": pKaWorkflow,
    "redox_potential": RedoxPotentialWorkflow,
    "scan": ScanWorkflow,
    "spin_states": SpinStatesWorkflow,
    "tautomers": TautomerWorkflow,
}
