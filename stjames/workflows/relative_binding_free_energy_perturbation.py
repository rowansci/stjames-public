"""Workflows covering RBFE graph construction and endpoint FEP execution."""

from typing import Annotated, Literal, Self

from pydantic import AfterValidator, PositiveInt, model_validator

from ..base import Base, round_float, round_optional_float
from ..message import Message
from ..method import Method
from ..molecule import Molecule
from ..pdb import PDB
from ..types import UUID
from .workflow import Workflow


class TMDRBFESettings(Base):
    """
    TMD-specific simulation parameters shared across RBFE graph edges.

    This settings bundle targets the TMD runner; other engines should define
    their own settings model if their controls differ.

    :param forcefield: Registered St. James method corresponding to a force field.
    :param n_eq_steps: Equilibration steps per lambda window.
    :param n_frames: Production frames saved per lambda window.
    :param steps_per_frame: MD integration steps per saved frame.
    :param n_windows: Maximum number of lambda windows considered for bisection.
    :param min_overlap: Minimum acceptable overlap during schedule bisection.
    :param target_overlap: Desired overlap after HREX optimization.
    :param water_sampling_padding: Extra nanometers added to the solvent sampling radius.
    :param rest_max_temperature_scale: Maximum effective temperature scaling for REST.
    :param rest_temperature_scale_interpolation: Functional form used for REST scaling.
    :param local_md_steps: Number of local MD steps per frame (`0` disables local MD).
    :param local_md_k: Spring constant used during local MD.
    :param local_md_radius: Sphere radius in nanometers for the local MD region.
    :param local_md_free_reference: Whether to free the reference frame during local MD.
    """

    forcefield: Method = Method.SMIRNOFF_2_2_1_AMBER_AM1BCC
    n_eq_steps: PositiveInt = 200_000
    n_frames: PositiveInt = 2_000
    steps_per_frame: PositiveInt = 400
    n_windows: PositiveInt = 48
    min_overlap: float = 0.667
    target_overlap: float = 0.667
    water_sampling_padding: float = 0.4
    rest_max_temperature_scale: float = 1.0
    rest_temperature_scale_interpolation: Literal["exponential", "linear"] = "exponential"
    local_md_steps: int = 390
    local_md_k: float = 10_000.0
    local_md_radius: float = 1.2
    local_md_free_reference: bool = False


class RBFEResult(Base):
    """
    Aggregate RBFE outcome for a single ligand.

    :param dg: Predicted binding free energy difference (kcal/mol).
    :param dg_err: Uncertainty estimate on `dg`.
    """

    dg: Annotated[float, AfterValidator(round_float(3))]
    dg_err: Annotated[float, AfterValidator(round_float(3))]


class RBFEGraphEdge(Base):
    """
    RBFE Edge definition with optional FEP edge results.

    :param mol_a: Source ligand identifier.
    :param mol_b: Target ligand identifier.
    :param core: Atom-mapping pairs describing the shared core, if any.
    :param score: Optional score used during graph construction.
    :param complex_dg: Predicted complex-leg free energy difference (kcal/mol).
    :param complex_dg_err: Uncertainty on `complex_dg`.
    :param solvent_dg: Predicted solvent-leg free energy difference (kcal/mol).
    :param solvent_dg_err: Uncertainty on `solvent_dg`.
    :param vacuum_dg: Predicted vacuum-leg free energy difference (kcal/mol).
    :param vacuum_dg_err: Uncertainty on `vacuum_dg`.
    :param ddg: Combined cycle result derived from complex and solvent legs.
    :param ddg_err: Uncertainty on `ddg`.
    """

    ligand_a: str
    ligand_b: str
    core: list[tuple[int, int]] | None = None
    score: float | None = None

    complex_dg: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    complex_dg_err: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    solvent_dg: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    solvent_dg_err: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    vacuum_dg: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    vacuum_dg_err: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    ddg: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    ddg_err: Annotated[float | None, AfterValidator(round_optional_float(3))] = None


class RBFEGraph(Base):
    """
    Minimal RBFE graph container.

    :param edges: Directed edges describing the ligand pairs to simulate.
    """

    edges: list[RBFEGraphEdge]


class RBFEGraphWorkflow(Workflow):
    """
    Workflow that builds an rbfe graph for a set of ligands.

    :param ligands: Mapping from ligand identifiers to `Molecule` objects.
    :param graph: Optional RBFE graph output populated after the build step.
    :param mode: Graph construction strategy (`"greedy"` or `"star_map"`).
    :param hub_compound_id: Ligand identifier to serve as the hub when `mode="star_map"`.
    :param refine_cutoff: Optional cutoff used to re-run atom mapping refinement.
    :param greedy_scoring: Edge scoring heuristic for greedy mode.
    :param greedy_k_min_cut: Target edge-connectivity (`k`) for greedy augmentation.
    """

    ligands: dict[str, Molecule]
    graph: RBFEGraph | None = None

    mode: Literal["greedy", "star_map"] = "greedy"
    hub_compound_id: str | None = None
    refine_cutoff: float | None = None
    greedy_scoring: Literal["best", "jaccard", "dummy_atoms"] = "best"
    greedy_k_min_cut: PositiveInt = 3

    @model_validator(mode="after")
    def validate_builder(self) -> Self:
        """
        Validate the builder inputs.

        :raises ValueError: If star-map mode omits `hub_compound_id` or fewer than two ligands are provided.
        """
        if self.mode == "star_map" and not self.hub_compound_id:
            raise ValueError("hub_compound_id is required when mode='star_map'")
        if len(self.ligands) < 2:
            raise ValueError("Provide at least two ligands to build an RBFE graph")
        return self


class RBFEDiagnostics(Base):
    """
    Quality-control metrics gathered during endpoint RBFE.

    :param cycle_closure_rms: RMS error across completed thermodynamic cycles.
    :param windows_completed: Count of successfully converged lambda windows.
    :param windows_failed: Count of failed lambda windows.
    :param notes: Structured messages describing noteworthy events.
    """

    cycle_closure_rms: Annotated[float | None, AfterValidator(round_optional_float(3))] = None
    windows_completed: PositiveInt | None = None
    windows_failed: PositiveInt | None = None
    notes: list[Message] = []


class RelativeBindingFreeEnergyPerturbationWorkflow(Workflow):
    """
    Workflow for running relative binding free energy perturbation simulations.

    :param ligands: Mapping from ligand identifiers to `Molecule` objects.
    :param graph: RBFE graph topology.
    :param target: PDB object or the UUID of the PDB object used for simulation.
    :param ligand_dg_results: Optional per-ligand FEP summaries produced downstream.
    :param diagnostics: Optional aggregate QC metrics.
    :param settings: Simulation controls shared across all RBFE edges.
    """

    ligands: dict[str, Molecule]
    graph: RBFEGraph
    target: PDB | UUID

    settings: TMDRBFESettings

    ligand_dg_results: dict[str, RBFEResult] | None = None
    diagnostics: RBFEDiagnostics | None = None
