from typing import Annotated, Any

from pydantic import AfterValidator, PositiveFloat, PositiveInt, model_validator

from ..base import Base, round_float
from ..pdb import PDB
from ..types import UUID, round_list
from .workflow import ProteinStructureWorkflow, SMILESWorkflow


class BindingPoseContact(Base):
    """
    A single protein–ligand contact from an MD trajectory.

    :param protein_atom_index: the index of the protein atom
    :param ligand_atom_index: the index of the ligand atom
    :occupancy: the probability of seeing this interaction in a frame, between 0 and 1
    """

    protein_atom_index: int
    ligand_atom_index: int
    occupancy: Annotated[float, AfterValidator(round_float(3))]


class BindingPoseTrajectory(Base):
    """
    Represents a single trajectory looking at a binding pose.

    :param uuid: the UUID of the trajectory
    :param ligand_rmsd: the RMSD of the ligand vs. starting pose (aligning the protein)
    :param contacts: the conserved binding-pose contacts
    """

    uuid: UUID
    ligand_rmsd: Annotated[list[float], AfterValidator(round_list(3))] = []
    contacts: list[BindingPoseContact] = []


class ProteinMDSettingsMixin(Base):
    """
    Mix-in for various settings used in running protein MD.

    :param equilibration_time_ns: how long to equilibrate trajectories for, in nanoseconds
    :param simulation_time_ns: how long to run trajectories for, in nanoseconds
    :param temperature: the temperature, in K
    :param pressure_atm: the pressure, in atm
    :param langevin_timescale_ps: the timescale for the Langevin integrator, in inverse picoseconds
    :param timestep_fs: the timestep, in femtoseconds
    :param constrain_hydrogens: whether or not to use SHAKE to freeze bonds to hydrogen
    :param nonbonded_cutoff: the nonbonded cutoff for particle-mesh Ewald, in Å
    :param protein_prune_cutoff: the cutoff past which residues will be deleted, in Å
    :param protein_restraint_cutoff: the cutoff past which alpha-carbons will be constrained, in Å
    :param protein_restraint_constant: the force constant for backbone restraints, in kcal/mol/Å**2
    :param ionic_strength_M: the ionic strength of the solution, in M (molar)
    :param water_buffer: the amount of water to add around the protein, in Å
    """

    equilibration_time_ns: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 1
    simulation_time_ns: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 10

    temperature: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 300
    pressure_atm: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 1.0
    langevin_timescale_ps: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 1.0

    timestep_fs: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 2
    constrain_hydrogens: bool = True
    nonbonded_cutoff: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 8.0

    protein_restraint_cutoff: Annotated[PositiveFloat, AfterValidator(round_float(3))] | None = None
    protein_restraint_constant: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 100

    ionic_strength_M: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 0.10
    water_buffer: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 6.0


class ProteinMolecularDynamicsWorkflow(ProteinMDSettingsMixin, ProteinStructureWorkflow):
    """
    Protein molecular dynamics workflow.

    Inherited:
    :param protein: PDB or UUID of the (holo) protein.
    :param equilibration_time_ns: how long to equilibrate trajectories for, in nanoseconds
    :param simulation_time_ns: how long to run trajectories for, in nanoseconds
    :param temperature: the temperature, in K
    :param pressure_atm: the pressure, in atm
    :param langevin_timescale_ps: the timescale for the Langevin integrator, in inverse picoseconds
    :param timestep_fs: the timestep, in femtoseconds
    :param constrain_hydrogens: whether or not to use SHAKE to freeze bonds to hydrogen
    :param nonbonded_cutoff: the nonbonded cutoff for particle-mesh Ewald, in Å
    :param protein_prune_cutoff: the cutoff past which residues will be deleted, in Å
    :param protein_restraint_cutoff: the cutoff past which alpha-carbons will be constrained, in Å
    :param protein_restraint_constant: the force constant for backbone restraints, in kcal/mol/Å**2
    :param ionic_strength_M: the ionic strength of the solution, in M (molar)
    :param water_buffer: the amount of water to add around the protein, in Å

    New:
    :param num_trajectories: the number of trajectories to run
    :param save_solvent: whether solvent should be saved

    Results:
    :param minimized_protein_uuid: UUID of final system PDB
    :param bonds: which atoms are bonded to which other atoms
    :param trajectories: the UUID for each trajectory
    """

    num_trajectories: PositiveInt = 1
    save_solvent: bool = False

    minimized_protein_uuid: UUID | None = None
    bonds: list[tuple[int, int]] = []
    trajectories: list[UUID] = []


class PoseAnalysisMolecularDynamicsWorkflow(ProteinMDSettingsMixin, ProteinStructureWorkflow, SMILESWorkflow):
    """
    Pose-analysis molecular dynamics workflow.

    Inherited:
    :param initial_smiles: ligand's SMILES
    :param protein: PDB or UUID of the (holo) protein.
    :param equilibration_time_ns: how long to equilibrate trajectories for, in nanoseconds
    :param simulation_time_ns: how long to run trajectories for, in nanoseconds
    :param temperature: the temperature, in K
    :param pressure_atm: the pressure, in atm
    :param langevin_timescale_ps: the timescale for the Langevin integrator, in inverse picoseconds
    :param timestep_fs: the timestep, in femtoseconds
    :param constrain_hydrogens: whether or not to use SHAKE to freeze bonds to hydrogen
    :param nonbonded_cutoff: the nonbonded cutoff for particle-mesh Ewald, in Å
    :param protein_prune_cutoff: the cutoff past which residues will be deleted, in Å
    :param protein_restraint_cutoff: the cutoff past which alpha-carbons will be constrained, in Å
    :param protein_restraint_constant: the force constant for backbone restraints, in kcal/mol/Å**2
    :param ionic_strength_M: the ionic strength of the solution, in M (molar)
    :param water_buffer: the amount of water to add around the protein, in Å

    New:
    :param protein_uuid: UUID of the (holo) protein. DEPRECATED.
    :param ligand_residue_name: ligand's residue name
    :param num_trajectories: the number of trajectories to run
    :param save_solvent: whether solvent should be saved

    Results:
    :param minimized_protein_uuid: UUID of final system PDB
    :param bonds: which atoms are bonded to which other atoms
    :param trajectories: for each replicate, a UUID and the corresponding analysis results
    """

    protein_uuid: UUID | None = None
    ligand_residue_name: str = "LIG"

    num_trajectories: PositiveInt = 1
    save_solvent: bool = False

    minimized_protein_uuid: UUID | None = None
    bonds: list[tuple[int, int]] = []
    trajectories: list[BindingPoseTrajectory] = []

    @model_validator(mode="before")
    def harmonize_protein_uuid(cls, data: Any) -> Any:
        """
        Syncs data between "protein_uuid" and "protein" field.
        """

        if data.get("protein_uuid", False):
            data["protein"] = data["protein_uuid"]

        return data
