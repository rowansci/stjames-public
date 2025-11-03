"""Conformer search workflow."""

from abc import ABC
from typing import Annotated, Literal, Self, Sequence, TypeVar

from pydantic import AfterValidator, BaseModel, Field, PositiveInt, field_validator, model_validator

from ..base import Base, LowercaseStrEnum
from ..constraint import Constraint
from ..method import Method, XTBMethod
from ..mode import Mode
from ..molecule import Molecule
from ..settings import Settings
from ..types import UUID, FloatPerAtom, round_float_per_atom
from .multistage_opt import MultiStageOptMixin
from .workflow import MoleculeWorkflow, SMILESWorkflow

_sentinel = object()

_T = TypeVar("_T")
_U = TypeVar("_U")


def check_sentinel(value: _T, default: _U) -> _T | _U:
    """Return value unless _sentinel, then return default."""
    return default if value is _sentinel else value


class ScreeningSettings(BaseModel):
    """
    Settings for determining unique and useful conformers.

    :param energy_threshold: maximum relative energy for screening
    :param rotational_constants_threshold: maximum difference in rotational constants for screening
    :param rmsd: Cartesian RMSD for screening
    :param max_confs: maximum number of conformers to keep
    """

    energy_threshold: float | None = None  # kcal/mol
    rotational_constants_threshold: float | None = 0.02
    rmsd: float | None = 0.25
    max_confs: int | None = None


class ConformerClusteringDescriptor(LowercaseStrEnum):
    """
    Potential descriptors to employ in conformer clustering.
    """

    SOLVENT_ACCESSIBLE_SURFACE_AREA = "solvent_accessible_surface_area"
    POLAR_SOLVENT_ACCESSIBLE_SURACE_AREA = "polar_solvent_accessible_surface_area"
    RADIUS_OF_GYRATION = "radius_of_gyration"
    PLANE_OF_BEST_FIT = "plane_of_best_fit"
    NORMALIZED_PRINCIPAL_MOMENT_RATIO_1 = "normalized_principal_moment_ratio_1"
    NORMALIZED_PRINCIPAL_MOMENT_RATIO_2 = "normalized_principal_moment_ratio_2"


class ConformerClusteringSettings(Base):
    """
    Settings for clustering conformers based on their three-dimensional properties.

    The properties used for clustering by default are:
    - Solvent-accessible surface area
    - Polar solvent-accessible surface area
    - Radius of gyration
    - Plane of best fit
    - Normalized principal moment ratios 1 and 2

    Rowan uses k-means clustering to identify representative conformers.
    This loosely follows Wilcken and co-workers (10.1007/s10822-020-00337-7).

    :param num_clusters: the number of clusters to include
    :param conformers_per_cluster: the number of compounds to pick from each cluster
    """

    descriptors: list[ConformerClusteringDescriptor] = [
        ConformerClusteringDescriptor.SOLVENT_ACCESSIBLE_SURFACE_AREA,
        ConformerClusteringDescriptor.POLAR_SOLVENT_ACCESSIBLE_SURACE_AREA,
        ConformerClusteringDescriptor.RADIUS_OF_GYRATION,
        ConformerClusteringDescriptor.PLANE_OF_BEST_FIT,
        ConformerClusteringDescriptor.NORMALIZED_PRINCIPAL_MOMENT_RATIO_1,
        ConformerClusteringDescriptor.NORMALIZED_PRINCIPAL_MOMENT_RATIO_2,
    ]

    num_clusters: PositiveInt = 5
    conformers_per_cluster: PositiveInt = 3


class ConformerGenSettings(BaseModel):
    """
    Conformer generation settings.

    Conformers are generated and an initial screening is performed to remove duplicates and high-energy conformers.

    :param mode: Mode for calculations
    :param conf_opt_method: method for the optimization
    :param screening: post-generation screening settings
    :param constraints: constraints for conformer generation
    :param nci: add a constraining potential for non-covalent interactions
    :param max_confs: maximum number of conformers to keep
    """

    mode: Mode = Mode.RAPID
    conf_opt_method: XTBMethod = Method.GFN_FF
    screening: ScreeningSettings | None = None
    constraints: Sequence[Constraint] = tuple()
    nci: bool = False
    max_confs: int | None = None

    def __str__(self) -> str:
        """Return a string representation of the ConformerGenSettings."""
        return repr(self)

    def __repr__(self) -> str:
        """Return a string representation of the ConformerGenSettings."""
        return f"<{type(self).__name__} {self.mode.name}>"


class ETKDGSettings(ConformerGenSettings):
    """
    Settings for ETKDG conformer generation.

    Inherited:
    :param mode: Mode for calculations
    :param conf_opt_method: method for the optimization
    :param screening: post-generation screening settings
    :param constraints: constraints for conformer generation
    :param nci: add a constraining potential for non-covalent interactions (not supported in ETKDG)
    :param max_confs: maximum number of conformers to keep

    New:
    :param num_initial_confs: number of initial conformers to generate
    :param num_confs_considered: number of conformers to consider for optimization
    :param num_confs_taken: number of final conformers to take
    :param max_mmff_energy: MMFF energy cutoff
    :param max_mmff_iterations: MMFF optimization iterations
    """

    num_initial_confs: int = 300
    num_confs_considered: int = 100
    max_mmff_iterations: int = 500
    max_mmff_energy: float | None = 30
    settings_type: Literal["etkdg"] = "etkdg"

    @field_validator("constraints")
    def check_constraints(cls, constraints: Sequence[Constraint]) -> Sequence[Constraint]:
        if constraints:
            raise ValueError("ETKDG does not support constraints")

        return tuple(constraints)

    @field_validator("nci")
    def check_nci(cls, nci: bool) -> Literal[False]:
        if nci:
            raise ValueError("ETKDG does not support NCI")

        return nci

    @model_validator(mode="after")
    def validate_and_build(self) -> Self:
        match self.mode:
            case Mode.MANUAL:
                pass
            case Mode.RECKLESS:
                self.num_initial_confs = 200
                self.num_confs_considered = 50
                self.max_confs = 20 if self.max_confs is None else self.max_confs
                self.max_mmff_energy = 20
            case Mode.RAPID:
                self.max_confs = 50 if self.max_confs is None else self.max_confs
                self.conf_opt_method = Method.GFN0_XTB
            case _:
                raise NotImplementedError(f"Unsupported mode: {self.mode}")

        return self


class iMTDSpeeds(LowercaseStrEnum):
    MEGAQUICK = "megaquick"
    SUPERQUICK = "superquick"
    QUICK = "quick"
    NORMAL = "normal"
    EXTENSIVE = "extensive"


class iMTDSettings(ConformerGenSettings, ABC):
    """
    Settings for iMTD style conformer generation.

    RECKLESS:
        - GFN-FF//MTD(GFN-FF)
        - Megaquick
            - No GC
            - No rotamer metadynamics
            - Energy window = 5.0
            - Run scaling factor = 0.5
            - 6 MTD runs
    RAPID:
        - GFN0//MTD(GFN-FF)
        - Superquick
            - No GC
            - No rotamer metadynamics
            - Energy window = 5.0
            - Run scaling factor = 0.5
            - 6 MTD runs
    CAREFUL:
        - GFN2//MTD(GFN-FF)
        - Quick
            - GC (for iMTD-GC)
            - Rotamer metadynamics (for iMTD-GC)
            - Energy window = 5.0
            - Run scaling factor = 0.5
            - 6 MTD runs
    METICULOUS:
        - GFN2//MTD(GFN2)
        - "Normal"
            - GC (for iMTD-GC)
            - Rotamer metadynamics (for iMTD-GC)
            - Energy window = 6.0
            - Run scaling factor = 1
            - 14 MTD runs (2 with extreme values)


    See https://github.com/crest-lab/crest/blob/5ca82feb2ec4df30a0129db957163c934f085952/src/choose_settings.f90#L202
    and https://github.com/crest-lab/crest/blob/5ca82feb2ec4df30a0129db957163c934f085952/src/confparse.f90#L825
    for how quick, superquick, and megaquick are defined.

    Additional notes:
        Extensive mode
        - GC
        - Rotamer metadynamics
        - Energy window = 8.0
        - Run scaling factor = 2
        - 14 MTD runs (2 with extreme values)

        --NCI may switch things to QUICK?

    Inherited:
    :param mode: Mode for calculations
    :param conf_opt_method: method for the optimization
    :param screening: post-generation screening settings (not used)
    :param constraints: constraints to add
    :param nci: add an ellipsoide potential around the input structure
    :param max_confs: maximum number of conformers to keep

    New:
    :param mtd_method: method for the metadynamics
    :param speed: speed of the calculations (CREST specific setting)
    :param reopt: re-optimize conformers (corrects for the lack of rotamer metadynamics and GC)
    :param free_energy_weights: calculate frequencies and re-weight based on free energies
    """

    mtd_method: XTBMethod = Method.GFN_FF
    mtd_runtype: str = "imtd-gc"

    speed: iMTDSpeeds = iMTDSpeeds.QUICK
    reopt: bool = _sentinel  # type: ignore [assignment]
    free_energy_weights: bool = False
    settings_type: Literal["imtd"] = "imtd"

    @model_validator(mode="after")
    def validate_and_build_imtdgc_settings(self) -> Self:
        match self.mode:
            case Mode.MANUAL:
                if self.reopt is _sentinel:
                    raise ValueError("Must specify reopt with MANUAL mode")
            case Mode.RECKLESS:  # GFN-FF//MTD(GFN-FF)
                self.max_confs = 20 if self.max_confs is None else self.max_confs
                self.speed = iMTDSpeeds.MEGAQUICK
                self.reopt = check_sentinel(self.reopt, True)
            case Mode.RAPID:  # GFN0//MTD(GFN-FF)
                self.max_confs = 50 if self.max_confs is None else self.max_confs
                self.speed = iMTDSpeeds.SUPERQUICK
                self.conf_opt_method = Method.GFN0_XTB
                self.reopt = check_sentinel(self.reopt, True)
            case Mode.CAREFUL:  # GFN2//MTD(GFN-FF)
                self.speed = iMTDSpeeds.QUICK
                self.conf_opt_method = Method.GFN2_XTB
                self.reopt = check_sentinel(self.reopt, False)
            case Mode.METICULOUS:  # GFN2//MTD(GFN2)
                self.speed = iMTDSpeeds.NORMAL
                self.mtd_method = Method.GFN2_XTB
                self.conf_opt_method = Method.GFN2_XTB
                self.reopt = check_sentinel(self.reopt, False)
            # case Mode.EXTREME: # GFN2//MTD(GFN2)
            #     self.mtd_method = Method.GFN2_XTB
            #     self.conf_opt_method = Method.GFN2_XTB
            #     self.speed = iMTDSpeeds.EXTENSIVE
            #     self.reopt = check_sentinel(self.reopt, False)
            case _:
                raise NotImplementedError(f"Unsupported mode: {self.mode}")

        return self


class iMTDGCSettings(iMTDSettings):
    run_type: str = "imtdgc"


class iMTDsMTDSettings(iMTDSettings):
    run_type: str = "imtd-smtd"


class LyrebirdSettings(ConformerGenSettings):
    """
    Settings for Lyrebird-based conformer generation.

    Inherited:
    :param mode: Mode for calculations
    :param conf_opt_method: method for the optimization
    :param screening: post-generation screening settings
    :param constraints: constraints for conformer generation (not supported)
    :param nci: add a constraining potential for non-covalent interactions (not supported)
    :param max_confs: maximum number of conformers to keep

    New:
    :param num_initial_confs: number of initial conformers to generate
    """

    num_initial_confs: int = 300
    settings_type: Literal["lyrebird"] = "lyrebird"

    @field_validator("constraints")
    def check_constraints(cls, constraints: Sequence[Constraint]) -> Sequence[Constraint]:
        if constraints:
            raise ValueError("Lyrebird does not support constraints")

        return tuple(constraints)

    @field_validator("nci")
    def check_nci(cls, nci: bool) -> Literal[False]:
        if nci:
            raise ValueError("Lyrebird does not support NCI")

        return nci


class MonteCarloMultipleMinimumSettings(ConformerGenSettings):
    """
    Settings for Monte-Carlo-multiple-minimum-based conformer generation.
    Default values recommended by Nick Casetti.

    Inherited:
    :param mode: Mode for calculations
    :param conf_opt_method: method for the optimization
    :param screening: post-generation screening settings
    :param constraints: constraints for conformer generation (not supported)
    :param nci: add a constraining potential for non-covalent interactions (not supported)
    :param max_confs: maximum number of conformers to keep

    New:
    :param num_monte_carlo_iterations: number of Monte Carlo iterations to run
    :param rmsd_threshold: the threshold to determine if MCMM output structures are identical
    :param energy_window: maximum energy window above the minimum-energy conformer above which to retain (kcal/mol)
    :param monte_carlo_settings: the way the actual energy will be computed for the Monte-Carlo steps
    """

    energy_settings: Settings = Settings(method=Method.AIMNET2_WB97MD3)

    num_monte_carlo_iterations: int = 100
    rmsd_threshold: float = 0.5
    energy_window: float = 200.0

    settings_type: Literal["monte_carlo_multiple_minimum"] = "monte_carlo_multiple_minimum"

    @field_validator("constraints")
    def check_constraints(cls, constraints: Sequence[Constraint]) -> Sequence[Constraint]:
        if constraints:
            raise ValueError("MCMM does not support constraints")

        return tuple(constraints)

    @field_validator("nci")
    def check_nci(cls, nci: bool) -> Literal[False]:
        if nci:
            raise ValueError("MCMM does not support NCI")

        return nci


ConformerGenSettingsUnion = Annotated[ETKDGSettings | iMTDSettings | LyrebirdSettings | MonteCarloMultipleMinimumSettings, Field(discriminator="settings_type")]


class ConformerGenMixin(BaseModel):
    """
    Mixin for workflows that need conformer generation.

    :param conf_gen_mode: Mode for calculations
    :param conf_gen_settings: settings for conformer generation
    :param constraints: constraints to add
    :param nci: add a constraining potential for non-covalent interactions
    :param max_confs: maximum number of conformers to keep
    :param clustering_settings: how to cluster the conformers (if at all)
    """

    conf_gen_mode: Mode = Mode.RAPID
    conf_gen_settings: ConformerGenSettingsUnion = _sentinel  # type: ignore [assignment]
    constraints: Sequence[Constraint] = tuple()
    nci: bool = False
    max_confs: int | None = None

    conformer_clustering_settings: ConformerClusteringSettings | None = None

    @model_validator(mode="after")
    def validate_and_build_conf_gen_settings(self) -> Self:
        """Validate and build the ConformerGenSettings."""
        if self.conf_gen_settings is not _sentinel:
            return self

        match self.conf_gen_mode:
            case Mode.MANUAL:
                if self.conf_gen_settings is _sentinel:
                    raise ValueError("Must specify conf_gen_settings with MANUAL mode")

            case Mode.RECKLESS | Mode.RAPID:
                # ETKDGSettings will error if constraints or nci are set
                self.conf_gen_settings = ETKDGSettings(mode=self.conf_gen_mode, constraints=self.constraints, nci=self.nci, max_confs=self.max_confs)
            case Mode.CAREFUL | Mode.METICULOUS:
                self.conf_gen_settings = iMTDSettings(mode=self.conf_gen_mode, constraints=self.constraints, nci=self.nci, max_confs=self.max_confs)

            case _:
                raise NotImplementedError(f"Unsupported mode: {self.conf_gen_mode}")

        return self


class ConformerSearchMixin(ConformerGenMixin, MultiStageOptMixin):
    """
    Mixin for workflows that need conformer search—a combination of conformer generation and optimization.

    Inherited (ConformerGenMixin):
    :param conf_gen_mode: Mode for conformer generation
    :param mso_mode: Mode for MultiStageOptSettings
    :param conf_gen_settings: settings for conformer generation
    :param nci: add a constraining potential for non-covalent interactions

    Inherited (MultiStageOptMixin):
    :param multistage_opt_settings: set by mso_mode unless mode=MANUAL (ignores additional settings if set)
    :param solvent: solvent to use
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param transition_state: whether this is a transition state
    :param frequencies: whether to calculate frequencies

    Inherited (Both):
    :param constraints: constraints to add (diamond inheritance, works as expected)
    """

    def __str__(self) -> str:
        """Return a string representation of the ConformerSearch workflow."""
        return repr(self)

    def __repr__(self) -> str:
        """Return a string representation of the ConformerSearch workflow."""
        return f"<{type(self).__name__} {self.conf_gen_mode.name} {self.mso_mode.name}>"

    @model_validator(mode="after")
    def remove_ts_constraints(self) -> Self:
        """
        Remove constraints from optimization if a TS.

        Also affects Manual Mode.
        """
        msos = self.multistage_opt_settings
        if msos.transition_state and msos.constraints:
            msos.constraints = []
            for opt_set in msos.optimization_settings:
                opt_set.opt_settings.constraints = []

        return self


class ConformerSearchWorkflow(ConformerSearchMixin, SMILESWorkflow, MoleculeWorkflow):
    """
    ConformerSearch Workflow.

    This workflow supports both SMILES and 3D molecular input. Some conformer generation settings
    support both methods; others (like CREST) require 3D information. Only one should be supplied.

    Inherited:
    :param initial_molecule: Molecule of interest
    :param initial_smiles: SMILES of the molecule of interest
    :param conf_gen_mode: Mode for calculations
    :param conf_gen_settings: settings for conformer generation
    :param mso_mode: Mode for MultiStageOptSettings
    :param multistage_opt_settings: set by mode unless mode=MANUAL (ignores additional settings if set)
    :param solvent: solvent to use
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param constraints: constraints to add
    :param transition_state: whether this is a transition state
    :param frequencies: whether to calculate frequencies

    Ignored:
    :param mode: Mode to use (not used)

    New:
    :param conformer_uuids: list of UUIDs of the Molecules generated
    :param energies: energies of the molecules
    """

    initial_smiles: str = ""
    initial_molecule: Molecule | None = None  # type: ignore [assignment]

    # Results
    conformer_uuids: list[list[UUID | None]] = Field(default_factory=list)
    energies: Annotated[FloatPerAtom, AfterValidator(round_float_per_atom(6))] = Field(default_factory=list)

    @model_validator(mode="after")
    def validate_mol_input(self) -> Self:
        """Ensure that only one of initial_molecule or initial_smiles is set."""

        if not (bool(self.initial_smiles) ^ bool(self.initial_molecule)):
            raise ValueError("Can only set one of initial_molecule and initial_smiles")

        if isinstance(self.conf_gen_settings, iMTDSettings) and (self.initial_molecule is None):
            raise ValueError("iMTDSettings requires `initial_molecule` to be set")

        return self
