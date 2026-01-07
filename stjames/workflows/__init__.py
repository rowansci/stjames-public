# ruff: noqa: F405

from typing import Literal

from .admet import *
from .basic_calculation import *
from .batch_docking import *
from .bde import *
from .conformer import *
from .conformer_search import *
from .descriptors import *
from .docking import *
from .double_ended_ts_search import DoubleEndedTSSearchWorkflow
from .electronic_properties import *
from .fukui import *
from .hydrogen_bond_basicity import *
from .ion_mobility import *
from .irc import *
from .macropka import *
from .membrane_permeability import *
from .molecular_dynamics import *
from .msa import *
from .multistage_opt import *
from .nmr import *
from .pka import *
from .pocket_detection import Pocket, PocketDetectionWorkflow
from .pose_analysis_md import *
from .protein_binder_design import *
from .protein_cofolding import *
from .redox_potential import *
from .relative_binding_free_energy_perturbation import (
    RBFEGraphWorkflow,
    RelativeBindingFreeEnergyPerturbationWorkflow,
)
from .scan import *
from .solubility import *
from .solvent_dependent_conformers import *
from .spin_states import *
from .strain import *
from .tautomer import *
from .workflow import *

WORKFLOW_NAME = Literal[
    "admet",
    "basic_calculation",
    "batch_docking",
    "bde",
    "membrane_permeability",
    "conformers",
    "conformer_search",
    "descriptors",
    "docking",
    "double_ended_ts_search",
    "electronic_properties",
    "relative_binding_free_energy_perturbation",
    "fukui",
    "hydrogen_bond_basicity",
    "ion_mobility",
    "irc",
    "macropka",
    "molecular_dynamics",
    "msa",
    "multistage_opt",
    "nmr",
    "pka",
    "pocket_detection",
    "pose_analysis_md",
    "protein_cofolding",
    "protein_binder_design",
    "rbfe_graph",
    "redox_potential",
    "scan",
    "solubility",
    "solvent_dependent_conformers",
    "spin_states",
    "strain",
    "tautomers",
]

WORKFLOW_MAPPING: dict[WORKFLOW_NAME, Workflow] = {
    "admet": ADMETWorkflow,  # type: ignore [dict-item]
    "basic_calculation": BasicCalculationWorkflow,  # type: ignore [dict-item]
    "batch_docking": BatchDockingWorkflow,  # type: ignore [dict-item]
    "bde": BDEWorkflow,  # type: ignore [dict-item]
    "membrane_permeability": MembranePermeabilityWorkflow,  # type: ignore [dict-item]
    "conformers": ConformerWorkflow,  # type: ignore [dict-item]
    "conformer_search": ConformerSearchWorkflow,  # type: ignore [dict-item]
    "descriptors": DescriptorsWorkflow,  # type: ignore [dict-item]
    "docking": DockingWorkflow,  # type: ignore [dict-item]
    "double_ended_ts_search": DoubleEndedTSSearchWorkflow,  # type: ignore [dict-item]
    "electronic_properties": ElectronicPropertiesWorkflow,  # type: ignore [dict-item]
    "relative_binding_free_energy_perturbation": RelativeBindingFreeEnergyPerturbationWorkflow,  # type: ignore [dict-item]
    "fukui": FukuiIndexWorkflow,  # type: ignore [dict-item]
    "hydrogen_bond_basicity": HydrogenBondBasicityWorkflow,  # type: ignore [dict-item]
    "ion_mobility": IonMobilityWorkflow,  # type: ignore [dict-item]
    "irc": IRCWorkflow,  # type: ignore [dict-item]
    "macropka": MacropKaWorkflow,  # type: ignore [dict-item]
    "molecular_dynamics": MolecularDynamicsWorkflow,  # type: ignore [dict-item]
    "multistage_opt": MultiStageOptWorkflow,  # type: ignore [dict-item]
    "msa": MSAWorkflow,  # type: ignore [dict-item]
    "nmr": NMRSpectroscopyWorkflow,  # type: ignore [dict-item]
    "pka": pKaWorkflow,  # type: ignore [dict-item]
    "pocket_detection": PocketDetectionWorkflow,  # type: ignore [dict-item]
    "pose_analysis_md": PoseAnalysisMolecularDynamicsWorkflow,  # type: ignore [dict-item]
    "protein_cofolding": ProteinCofoldingWorkflow,  # type: ignore [dict-item]
    "protein_binder_design": ProteinBinderDesignWorkflow,  # type: ignore [dict-item]
    "rbfe_graph": RBFEGraphWorkflow,  # type: ignore [dict-item]
    "redox_potential": RedoxPotentialWorkflow,  # type: ignore [dict-item]
    "scan": ScanWorkflow,  # type: ignore [dict-item]
    "solubility": SolubilityWorkflow,  # type: ignore [dict-item]
    "solvent_dependent_conformers": SolventDependentConformersWorkflow,  # type: ignore [dict-item]
    "spin_states": SpinStatesWorkflow,  # type: ignore [dict-item]
    "strain": StrainWorkflow,  # type: ignore [dict-item]
    "tautomers": TautomerWorkflow,  # type: ignore [dict-item]
}
