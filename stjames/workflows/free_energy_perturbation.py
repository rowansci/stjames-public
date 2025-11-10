"""Workflows covering RBFE graph construction and endpoint FEP execution."""

from __future__ import annotations

from typing import Literal

from pydantic import Field, PositiveInt, field_validator, model_validator

from ..base import Base
from ..message import Message
from ..molecule import Molecule
from ..pdb import PDB
from .workflow import Workflow


class RBFEGraphEdge(Base):
    """
    Edge definition plus optional FEP result.

    :param mol_a: Source ligand identifier.
    :param mol_b: Target ligand identifier.
    :param core: Atom-mapping pairs describing the shared core, if any.
    :param score: Optional score used during graph construction.
    :param complex_dg: Predicted complex-leg free energy difference (kcal/mol).
    :param complex_dg_err: Uncertainty on ``complex_dg``.
    :param solvent_dg: Predicted solvent-leg free energy difference (kcal/mol).
    :param solvent_dg_err: Uncertainty on ``solvent_dg``.
    :param vacuum_dg: Predicted vacuum-leg free energy difference (kcal/mol).
    :param vacuum_dg_err: Uncertainty on ``vacuum_dg``.
    :param ddg: Combined cycle result derived from complex and solvent legs.
    :param ddg_err: Uncertainty on ``ddg``.
    """

    mol_a: str
    mol_b: str
    core: list[tuple[int, int]] | None = None
    score: float | None = None

    complex_dg: float | None = None
    complex_dg_err: float | None = None
    solvent_dg: float | None = None
    solvent_dg_err: float | None = None
    vacuum_dg: float | None = None
    vacuum_dg_err: float | None = None
    ddg: float | None = None
    ddg_err: float | None = None


class RBFEGraph(Base):
    """
    Minimal RBFE graph container shared between build and FEP workflows.

    :param edges: Directed edges describing the ligand pairs to simulate.
    """

    edges: list[RBFEGraphEdge]


class BuildRBFEGraphWorkflow(Workflow):
    """
    Workflow that mirrors the ``examples/build_rbfe_graph.py`` CLI.

    :param ligands: Mapping from ligand identifiers to ``Molecule`` objects.
    :param graph: Optional RBFE graph output populated after the build step.
    :param mode: Graph construction strategy (``"greedy"`` or ``"star_map"``).
    :param hub_compound_id: Ligand identifier to serve as the hub when ``mode="star_map"``.
    :param refine_cutoff: Optional cutoff used to re-run atom mapping refinement.
    :param greedy_scoring: Edge scoring heuristic for greedy mode.
    :param greedy_k_min_cut: Target edge-connectivity (``k``) for greedy augmentation.
    :param report_interval: How often to report progress when iterating ligand pairs.
    :param verbose: Whether to emit extra diagnostics during graph generation.
    """

    ligands: dict[str, Molecule]
    graph: RBFEGraph | None = None

    mode: Literal["greedy", "star_map"] = "greedy"
    hub_compound_id: str | None = None
    refine_cutoff: float | None = None
    greedy_scoring: Literal["best", "jaccard", "dummy_atoms"] = "best"
    greedy_k_min_cut: PositiveInt = 3
    report_interval: PositiveInt = 100
    verbose: bool = False

    @model_validator(mode="after")
    def validate_builder(self) -> "BuildRBFEGraphWorkflow":
        """
        Validate the builder inputs.

        :raises ValueError: If star-map mode omits ``hub_compound_id`` or fewer than two ligands are provided.
        """
        if self.mode == "star_map" and not self.hub_compound_id:
            raise ValueError("hub_compound_id is required when mode='star_map'")
        if len(self.ligands) < 2:
            raise ValueError("Provide at least two ligands to build an RBFE graph")
        return self


class FEPDiagnostics(Base):
    """
    Quality-control metrics gathered during endpoint FEP.

    :param cycle_closure_rms_kcal_mol: RMS error across completed thermodynamic cycles.
    :param windows_completed: Count of successfully converged lambda windows.
    :param windows_failed: Count of failed lambda windows.
    :param notes: Structured messages describing noteworthy events.
    """

    cycle_closure_rms_kcal_mol: float | None = None
    windows_completed: PositiveInt | None = None
    windows_failed: PositiveInt | None = None
    notes: list[Message] = []


class FreeEnergyPerturbationWorkflow(Workflow):
    """
    Workflow that mirrors the ``examples/run_rbfe_graph.py`` CLI.

    :param ligands: Mapping from ligand identifiers to ``Molecule`` objects.
    :param graph: RBFE graph topology (typically the output of ``BuildRBFEGraphWorkflow``).
    :param pdb_structure: Prepared complex structure required for complex-leg simulations.

    :param ligand_dg_results: Optional per-ligand FEP summaries produced downstream.
    :param diagnostics: Optional aggregate QC metrics.

    :param forcefield: Serialized force field identifier loaded by the runner.
    :param legs: Ordered list of legs (vacuum/solvent/complex) to execute.
    :param n_eq_steps: Equilibration steps per lambda window.
    :param n_frames: Production frames saved per lambda window.
    :param steps_per_frame: MD integration steps per saved frame.
    :param n_windows: Maximum number of lambda windows considered for bisection.
    :param min_overlap: Minimum acceptable overlap during schedule bisection.
    :param target_overlap: Desired overlap after HREX optimization.
    :param water_sampling_padding: Extra nanometers added to the solvent sampling radius.
    :param rest_max_temperature_scale: Maximum effective temperature scaling for REST.
    :param rest_temperature_scale_interpolation: Functional form used for REST scaling.
    :param local_md_steps: Number of local MD steps per frame (``0`` disables local MD).
    :param local_md_k: Spring constant used during local MD.
    :param local_md_radius: Sphere radius in nanometers for the local MD region.
    :param local_md_free_reference: Whether to free the reference frame during local MD.
    :param store_trajectories: Whether to persist per-leg lambda endpoint trajectories.
    :param force_overwrite: Whether to overwrite an existing leg directory.
    :param experimental_field: Molecule property containing experimental measurements.
    :param experimental_units: Units corresponding to ``experimental_field``.
    :param seed: RNG seed applied to ligand shuffling and MD initialization.
    :param mps_workers: Number of MPS processes launched per GPU.
    :param n_gpus: Number of GPUs allocated for the workload (``None`` ⇒ auto-detect).
    :param output_dir: Directory to which CLI outputs should be written, if desired.
    """

    ligands: dict[str, Molecule]
    graph: RBFEGraph
    pdb_structure: PDB

    ligand_dg_results: dict[str, tuple[float, float]]
    diagnostics: FEPDiagnostics | None = None

    forcefield: str = "smirnoff_2_2_1_amber_am1bcc"
    legs: list[str] = Field(default_factory=lambda: ["vacuum", "solvent", "complex"])
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
    store_trajectories: bool = False
    force_overwrite: bool = False
    experimental_field: str = "kcal/mol experimental dG"
    experimental_units: Literal["kcal/mol", "kJ/mol", "uM", "nM"] = "kcal/mol"
    seed: int = 2025
    mps_workers: PositiveInt = 1
    n_gpus: PositiveInt | None = None
    output_dir: str | None = None

    @field_validator("legs")
    @classmethod
    def validate_legs(cls, legs: list[str]) -> list[str]:
        """
        Ensure the requested legs are valid and non-empty.

        :param legs: Ordered list of legs supplied to the workflow.
        :raises ValueError: If the list is empty or contains unsupported leg names.
        :return: The validated list of legs.
        """
        if not legs:
            raise ValueError("At least one leg must be specified")
        invalid = sorted(set(legs) - {"vacuum", "solvent", "complex"})
        if invalid:
            raise ValueError(f"Unsupported leg(s): {', '.join(invalid)}")
        return legs

    @model_validator(mode="after")
    def validate_inputs(self) -> "FreeEnergyPerturbationWorkflow":
        """
        Validate ligand, graph, and structure inputs before running FEP.

        :raises ValueError: If ligands or graph data are missing, empty, or incompatible with the requested legs.
        :return: The validated workflow instance.
        """
        if not self.ligands:
            raise ValueError("Provide at least one ligand")
        if self.graph is None:
            raise ValueError("graph must be provided")
        if not self.graph.edges:
            raise ValueError("graph must contain at least one edge")
        if "complex" in self.legs and not self.pdb_structure:
            raise ValueError("pdb_structure is required when running the complex leg")
        return self
