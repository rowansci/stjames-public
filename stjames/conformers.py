from abc import ABC
from typing import Annotated, Literal, Self, Sequence

from pydantic import AfterValidator, BaseModel, Field, NonNegativeFloat, PositiveFloat, PositiveInt, field_validator, model_validator

from .base import LowercaseStrEnum, round_float
from .constraint import Constraint
from .method import Method, XTBMethod
from .mode import Mode
from .settings import Settings
from .solvent import SolventModel, SolventSettings
from .workflows.multistage_opt import MultiStageOptMixin


class ConformerProperties(BaseModel):
    """
    Descriptors of a conformer's properties.

    :param solvent_accessible_surface_area: average SASA (Å²)
    :param polar_solvent_accessible_surface_area: average SASA for non-C/H elements (Å²)
    :param radius_of_gyration: radius of gyration (Å)
    """

    solvent_accessible_surface_area: Annotated[PositiveFloat, AfterValidator(round_float(3))]
    polar_solvent_accessible_surface_area: Annotated[NonNegativeFloat, AfterValidator(round_float(3))]
    radius_of_gyration: Annotated[PositiveFloat, AfterValidator(round_float(3))]


class ConformerClusteringDescriptor(LowercaseStrEnum):
    """Potential descriptors to employ in conformer clustering."""

    SOLVENT_ACCESSIBLE_SURFACE_AREA = "solvent_accessible_surface_area"
    POLAR_SOLVENT_ACCESSIBLE_SURACE_AREA = "polar_solvent_accessible_surface_area"
    RADIUS_OF_GYRATION = "radius_of_gyration"
    PLANE_OF_BEST_FIT = "plane_of_best_fit"
    NORMALIZED_PRINCIPAL_MOMENT_RATIO_1 = "normalized_principal_moment_ratio_1"
    NORMALIZED_PRINCIPAL_MOMENT_RATIO_2 = "normalized_principal_moment_ratio_2"


class ConformerClusteringSettings(BaseModel):
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

    :param num_clusters: number of clusters to include
    :param conformers_per_cluster: number of compounds to pick from each cluster
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

    :param screening: post-generation screening settings
    :param constraints: constraints for conformer generation
    :param nci: add a constraining potential for non-covalent interactions
    :param max_confs: maximum number of conformers to keep
    """

    constraints: Sequence[Constraint] = ()
    nci: bool = False
    max_confs: PositiveInt | None = None

    def __str__(self) -> str:
        """Return a string representation of the ConformerGenSettings."""
        return repr(self)

    def __repr__(self) -> str:
        """Return a string representation of the ConformerGenSettings."""
        extra = ""
        if self.constraints:
            extra += f" constraints={len(self.constraints)}"
        if self.nci:
            extra += " nci "
        if self.max_confs:
            extra += f" max_confs={self.max_confs}"

        return f"<{type(self).__name__} {extra}>"


class ETKDGSettings(ConformerGenSettings):
    """
    Settings for ETKDG conformer generation.

    Inherited:
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

    @classmethod
    def from_mode(cls, mode: Mode) -> Self:
        match mode:
            case Mode.RECKLESS:
                num_initial_confs = 200
                num_confs_considered = 50
                max_mmff_iterations = 500
                max_mmff_energy = 30
                max_confs = 20
                max_mmff_energy = 20
            case Mode.RAPID:
                num_initial_confs = 300
                num_confs_considered = 100
                max_mmff_iterations = 500
                max_mmff_energy = 30
                max_confs = 50
                max_mmff_energy = 30
            case _:
                raise NotImplementedError(f"Unsupported mode: {mode}")

        return cls(
            num_initial_confs=num_initial_confs,
            num_confs_considered=num_confs_considered,
            max_mmff_iterations=max_mmff_iterations,
            max_mmff_energy=max_mmff_energy,
            max_confs=max_confs,
        )


class iMTDSpeeds(LowercaseStrEnum):
    MEGAQUICK = "megaquick"
    SUPERQUICK = "superquick"
    QUICK = "quick"
    NORMAL = "normal"
    EXTENSIVE = "extensive"


class iMTDSettings(ConformerGenSettings, ABC):
    """
    Settings for iMTD style conformer generation.

    See https://github.com/crest-lab/crest/blob/5ca82feb2ec4df30a0129db957163c934f085952/src/choose_settings.f90#L202
    and https://github.com/crest-lab/crest/blob/5ca82feb2ec4df30a0129db957163c934f085952/src/confparse.f90#L825
    for how quick, superquick, and megaquick are defined.

    See build_imtd_setings(mode) for sensible defaults.

    Inherited:
    :param screening: post-generation screening settings (not used)
    :param constraints: constraints to add
    :param nci: add an ellipsoide potential around the input structure
    :param max_confs: maximum number of conformers to keep

    New:
    :param mtd_method: method for the metadynamics
    :param mtd_runtype: algorithm used
    :param speed: speed of the calculations (CREST specific setting)
    :param reopt: re-optimize conformers (corrects for the lack of rotamer metadynamics and GC)
    :param free_energy_weights: calculate frequencies and re-weight based on free energies
    :param energy_window: energy window used, in kcal/mol (CREST specific setting). if set, overrides default from speed
    :param solvent_settings: solvent to use, if any
    """

    settings_type: Literal["imtd"] = "imtd"

    mtd_method: XTBMethod = Method.GFN_FF
    mtd_runtype: str = "imtd-gc"
    speed: iMTDSpeeds = iMTDSpeeds.QUICK
    reopt: bool = False
    free_energy_weights: bool = False
    energy_window: Annotated[PositiveFloat, AfterValidator(round_float(1))] | None = None
    solvent_settings: SolventSettings | None = None

    @model_validator(mode="after")
    def validate_and_build_imtdgc_settings(self) -> Self:
        if self.solvent_settings and self.solvent_settings.model not in {SolventModel.ALPB, SolventModel.GBSA}:
            raise ValueError("Only ALPB or GBSA solvent models supported for iMTD conformer search!")

        return self

    @classmethod
    def from_mode(cls, mode: Mode) -> Self:
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
        """
        match mode:
            case Mode.RECKLESS:  # GFN-FF//MTD(GFN-FF)
                mtd_method = Method.GFN_FF
                speed = iMTDSpeeds.MEGAQUICK
                reopt = True
                max_confs: int | None = 20
            case Mode.RAPID:  # GFN0//MTD(GFN-FF)
                mtd_method = Method.GFN_FF
                speed = iMTDSpeeds.SUPERQUICK
                reopt = True
                max_confs = 50
            case Mode.CAREFUL:  # GFN2//MTD(GFN-FF)
                mtd_method = Method.GFN_FF
                speed = iMTDSpeeds.QUICK
                reopt = False
                max_confs = None
            case Mode.METICULOUS:  # GFN2//MTD(GFN2)
                mtd_method = Method.GFN2_XTB
                speed = iMTDSpeeds.NORMAL
                reopt = False
                max_confs = None
            case _:
                raise NotImplementedError(f"Unsupported mode: {mode}")

        return cls(
            mtd_method=mtd_method,
            speed=speed,
            reopt=reopt,
            max_confs=max_confs,
        )


class iMTDGCSettings(iMTDSettings):
    run_type: str = "imtdgc"


class iMTDsMTDSettings(iMTDSettings):
    run_type: str = "imtd-smtd"


class LyrebirdSettings(ConformerGenSettings):
    """
    Settings for Lyrebird-based conformer generation.

    Inherited:
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
    :param screening: post-generation screening settings
    :param constraints: constraints for conformer generation (not supported)
    :param nci: add a constraining potential for non-covalent interactions (not supported)
    :param max_confs: maximum number of conformers to keep

    New:
    :param num_monte_carlo_iterations: number of Monte Carlo iterations to run
    :param rmsd_threshold: threshold to determine if MCMM output structures are identical
    :param energy_window: maximum energy window above the minimum-energy conformer above which to retain (kcal/mol)
    :param monte_carlo_settings: energy computation method for Monte-Carlo steps
    """

    energy_settings: Settings = Settings(method=Method.AIMNET2_WB97MD3)

    num_monte_carlo_iterations: int = 250
    rmsd_threshold: float = 0.5
    energy_window: float = 20

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
    Mixin for classes need conformer generation.

    :param conf_gen_settings: settings for conformer generation
    :param constraints: constraints to add
    :param nci: add a constraining potential for non-covalent interactions
    :param max_confs: maximum number of conformers to keep
    :param clustering_settings: how to cluster the conformers (if at all)
    """

    conf_gen_settings: ConformerGenSettingsUnion | None
    constraints: Sequence[Constraint] = ()
    nci: bool = False
    max_confs: int | None = None

    conformer_clustering_settings: ConformerClusteringSettings | None = None


class ConformerSearchMixin(ConformerGenMixin, MultiStageOptMixin):
    """
    Mixin for classes need conformer search—a combination of conformer generation and optimization.

    Inherited (ConformerGenMixin):
    :param conf_gen_settings: settings for conformer generation
    :param nci: add a constraining potential for non-covalent interactions

    Inherited (MultiStageOptMixin):
    :param multistage_opt_settings: settings for the optimization
    :param solvent: solvent to use
    :param xtb_preopt: pre-optimize with xtb
    :param transition_state: whether this is a transition state
    :param frequencies: whether to calculate frequencies

    Inherited (Both):
    :param constraints: constraints to add (diamond inheritance, works as expected)
    """

    def __str__(self) -> str:
        """Return a string representation of the ConformerSearch workflow."""
        return repr(self)

    @model_validator(mode="after")
    def remove_ts_constraints(self) -> Self:
        """Remove constraints from optimization if a TS."""
        msos = self.multistage_opt_settings
        if msos.transition_state and msos.constraints:
            msos.constraints = []
            for opt_set in msos.optimization_settings:
                opt_set.opt_settings.constraints = []

        return self
