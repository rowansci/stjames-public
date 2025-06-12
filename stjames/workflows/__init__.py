# ruff: noqa: F405

from typing import Literal

from .admet import *
from .basic_calculation import *
from .bde import *
from .conformer import *
from .conformer_search import *
from .descriptors import *
from .docking import *
from .electronic_properties import *
from .fukui import *
from .hydrogen_bond_basicity import *
from .ion_mobility import *
from .irc import *
from .macropka import *
from .molecular_dynamics import *
from .multistage_opt import *
from .pka import *
from .pose_analysis_md import *
from .protein_cofolding import *
from .redox_potential import *
from .scan import *
from .solubility import *
from .spin_states import *
from .tautomer import *
from .workflow import *

WORKFLOW_NAME = Literal[
    "admet",
    "basic_calculation",
    "bde",
    "conformers",
    "conformer_search",
    "descriptors",
    "docking",
    "electronic_properties",
    "fukui",
    "hydrogen_bond_basicity",
    "ion_mobility",
    "irc",
    "macropka",
    "molecular_dynamics",
    "multistage_opt",
    "pka",
    "pose_analysis_md",
    "protein_cofolding",
    "redox_potential",
    "scan",
    "solubility",
    "spin_states",
    "tautomers",
]

WORKFLOW_MAPPING: dict[WORKFLOW_NAME, Workflow] = {
    "admet": ADMETWorkflow,  # type: ignore [dict-item]
    "basic_calculation": BasicCalculationWorkflow,  # type: ignore [dict-item]
    "bde": BDEWorkflow,  # type: ignore [dict-item]
    "conformers": ConformerWorkflow,  # type: ignore [dict-item]
    "conformer_search": ConformerSearchWorkflow,  # type: ignore [dict-item]
    "descriptors": DescriptorsWorkflow,  # type: ignore [dict-item]
    "docking": DockingWorkflow,  # type: ignore [dict-item]
    "electronic_properties": ElectronicPropertiesWorkflow,  # type: ignore [dict-item]
    "fukui": FukuiIndexWorkflow,  # type: ignore [dict-item]
    "hydrogen_bond_basicity": HydrogenBondBasicityWorkflow,  # type: ignore [dict-item]
    "ion_mobility": IonMobilityWorkflow,  # type: ignore [dict-item]
    "irc": IRCWorkflow,  # type: ignore [dict-item]
    "macropka": MacropKaWorkflow,  # type: ignore [dict-item]
    "molecular_dynamics": MolecularDynamicsWorkflow,  # type: ignore [dict-item]
    "multistage_opt": MultiStageOptWorkflow,  # type: ignore [dict-item]
    "pka": pKaWorkflow,  # type: ignore [dict-item]
    "pose_analysis_md": PoseAnalysisMolecularDynamicsWorkflow,  # type: ignore [dict-item]
    "protein_cofolding": ProteinCofoldingWorkflow,  # type: ignore [dict-item]
    "redox_potential": RedoxPotentialWorkflow,  # type: ignore [dict-item]
    "scan": ScanWorkflow,  # type: ignore [dict-item]
    "solubility": SolubilityWorkflow,  # type: ignore [dict-item]
    "spin_states": SpinStatesWorkflow,  # type: ignore [dict-item]
    "tautomers": TautomerWorkflow,  # type: ignore [dict-item]
}
