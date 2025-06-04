from typing import Self

from pydantic import PositiveFloat, PositiveInt, model_validator

from ..pdb import PDB
from ..types import UUID
from .workflow import MoleculeWorkflow


class PoseAnalysisMolecularDynamicsWorkflow(MoleculeWorkflow):
    """
    Pose analysis molecular dynamics workflow.

    Note that the protein can be supplied either by UUID or raw PDB object.
    We anticipate that the former will dominate deployed usage, but the latter is handy for isolated testing.
    If, for whatever reason, the workflow is initialized with both a `target_uuid` and a `target`, the UUID will be ignored.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow (currently unused)

    New:
    :param protein: PDB of the protein.
    :param protein_uuid: UUID of the protein.
    :param num_trajectories: how many trajectories to run
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

    Results:
    :param trajectories: UUIDs corresponding to Trajectory objects
    :param minimized_protein_uuid: UUID of final system PDB
    """

    protein: PDB | None = None
    protein_uuid: UUID | None = None

    num_trajectories: PositiveInt = 4
    equilibration_time_ns: PositiveFloat = 5
    simulation_time_ns: PositiveFloat = 10

    temperature: PositiveFloat = 300
    pressure_atm: PositiveFloat = 1.0
    langevin_timescale_ps: PositiveFloat = 1.0

    timestep_fs: PositiveFloat = 2
    constrain_hydrogens: bool = True
    nonbonded_cutoff: PositiveFloat = 8.0

    protein_prune_cutoff: PositiveFloat = 9.0
    protein_restraint_cutoff: PositiveFloat = 7.0
    protein_restraint_constant: PositiveFloat = 100

    ionic_strength_M: PositiveFloat = 0.10
    water_buffer: PositiveFloat = 6.0

    minimized_protein_uuid: UUID | None = None
    trajectories: list[UUID] = []

    @model_validator(mode="after")
    def check_cutoff_sanity(self) -> Self:
        """Check if protein is provided."""
        if self.protein_prune_cutoff < self.protein_restraint_cutoff:
            raise ValueError("Pruning cutoff must be larger than restraint cutoff")
        return self
