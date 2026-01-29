from typing import Annotated, Any

from pydantic import AfterValidator, PositiveFloat, PositiveInt, model_validator

from ..base import Base, round_float
from ..types import UUID, round_list
from .workflow import ProteinStructureWorkflow, SMILESWorkflow


class BindingPoseContact(Base):
    """
    A single protein–ligand contact from an MD trajectory.

    :param protein_atom_index: index of protein atom
    :param ligand_atom_index: index of ligand atom
    :occupancy: the probability of seeing this interaction in a frame, between 0 and 1
    """

    protein_atom_index: int
    ligand_atom_index: int
    occupancy: Annotated[float, AfterValidator(round_float(3))]


class ProteinMDTrajectory(Base):
    """
    Represents a single protein MD trajectory.

    :param uuid: UUID of trajectory
    """

    uuid: UUID


class BindingPoseTrajectory(ProteinMDTrajectory):
    """
    Represents a single trajectory looking at a binding pose.

    Inherited:
    :param uuid: UUID of trajectory

    New:
    :param ligand_rmsd: RMSD of ligand vs. starting pose (aligning protein)
    :param contacts: conserved binding-pose contacts
    """

    ligand_rmsd: Annotated[list[float], AfterValidator(round_list(3))] = []
    contacts: list[BindingPoseContact] = []


class ProteinMDSettingsMixin(Base):
    """
    Mix-in for various settings used in running protein MD.

    :param equilibration_time_ns: how long to equilibrate trajectories for, in ns
    :param simulation_time_ns: how long to run trajectories for, in ns
    :param temperature: temperature, in K
    :param pressure_atm: pressure, in atm
    :param langevin_timescale_ps: timescale for the Langevin integrator, in ps⁻¹
    :param timestep_fs: timestep, in femtoseconds
    :param constrain_hydrogens: whether or not to use SHAKE to freeze bonds to hydrogen
    :param nonbonded_cutoff: nonbonded cutoff for particle-mesh Ewald, in Å
    :param ionic_strength_M: ionic strength of the solution, in M (molar)
    :param water_buffer: amount of water to add around the protein, in Å
    """

    equilibration_time_ns: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 1
    simulation_time_ns: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 10

    temperature: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 300
    pressure_atm: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 1.0
    langevin_timescale_ps: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 1.0

    timestep_fs: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 2
    constrain_hydrogens: bool = True
    nonbonded_cutoff: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 8.0

    ionic_strength_M: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 0.10
    water_buffer: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 6.0


class ProteinMolecularDynamicsWorkflow(ProteinMDSettingsMixin, ProteinStructureWorkflow):
    """
    Protein molecular dynamics workflow.

    Inherited:
    :param protein: PDB or UUID of the (holo) protein.
    :param equilibration_time_ns: how long to equilibrate trajectories for, in ns
    :param simulation_time_ns: how long to run trajectories for, in ns
    :param temperature: temperature, in K
    :param pressure_atm: pressure, in atm
    :param langevin_timescale_ps: timescale for Langevin integrator, in ps⁻¹
    :param timestep_fs: timestep, in femtoseconds
    :param constrain_hydrogens: whether to use SHAKE to freeze bonds to hydrogen
    :param nonbonded_cutoff: nonbonded cutoff for particle-mesh Ewald, in Å
    :param protein_restraint_cutoff: cutoff past which alpha-carbons constrained, in Å
    :param protein_restraint_constant: force constant for backbone restraints, in kcal/mol/Å²
    :param ionic_strength_M: ionic strength of solution, in M (molar)
    :param water_buffer: amount of water to add around protein, in Å

    New:
    :param num_trajectories: number of trajectories to run
    :param save_solvent: whether solvent should be saved

    Results:
    :param minimized_protein_uuid: UUID of final system PDB
    :param bonds: which atoms are bonded to which other atoms
    :param trajectories: UUID for each trajectory
    """

    num_trajectories: PositiveInt = 1
    save_solvent: bool = False

    minimized_protein_uuid: UUID | None = None
    bonds: list[tuple[int, int]] = []
    trajectories: list[ProteinMDTrajectory] = []


class PoseAnalysisMolecularDynamicsWorkflow(ProteinMDSettingsMixin, ProteinStructureWorkflow, SMILESWorkflow):
    """
    Pose-analysis molecular dynamics workflow.

    Inherited:
    :param initial_smiles: ligand's SMILES
    :param protein: PDB or UUID of the (holo) protein.
    :param equilibration_time_ns: how long to equilibrate trajectories for, in ns
    :param simulation_time_ns: how long to run trajectories for, in ns
    :param temperature: temperature, in K
    :param pressure_atm: pressure, in atm
    :param langevin_timescale_ps: timescale for the Langevin integrator, in ps⁻¹
    :param timestep_fs: timestep, in femtoseconds
    :param constrain_hydrogens: whether or not to use SHAKE to freeze bonds to hydrogen
    :param nonbonded_cutoff: nonbonded cutoff for particle-mesh Ewald, in Å
    :param ionic_strength_M: ionic strength of the solution, in M (molar)
    :param water_buffer: amount of water to add around the protein, in Å

    New:
    :param protein_uuid: UUID of the (holo) protein. DEPRECATED.
    :param ligand_residue_name: ligand's residue name
    :param num_trajectories: number of trajectories to run
    :param save_solvent: whether solvent should be saved
    :param protein_restraint_cutoff: cutoff past which alpha-carbons will be constrained, in Å
    :param protein_restraint_constant: force constant for backbone restraints, in kcal/mol/Å²

    Results:
    :param minimized_protein_uuid: UUID of final system PDB
    :param bonds: which atoms are bonded to which other atoms
    :param trajectories: for each replicate, a UUID and the corresponding analysis results
    """

    protein_uuid: UUID | None = None
    ligand_residue_name: str = "LIG"

    num_trajectories: PositiveInt = 1
    save_solvent: bool = False

    protein_restraint_cutoff: Annotated[PositiveFloat, AfterValidator(round_float(3))] | None = None
    protein_restraint_constant: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 100

    minimized_protein_uuid: UUID | None = None
    bonds: list[tuple[int, int]] = []
    trajectories: list[BindingPoseTrajectory] = []

    @model_validator(mode="before")
    def harmonize_protein_uuid(cls, data: Any) -> Any:  # noqa: N805
        """
        Syncs data between "protein_uuid" and "protein" field.
        """
        if "protein_uuid" in data:
            data["protein"] = data["protein_uuid"]

        return data
