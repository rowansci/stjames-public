"""Conformer search workflow."""

from abc import ABC
from typing import Annotated, Self, Sequence, TypeVar

from pydantic import AfterValidator, BaseModel, Field, field_validator, model_validator

from ..base import LowercaseStrEnum
from ..constraint import Constraint
from ..method import Method, XTBMethod
from ..mode import Mode
from ..types import UUID, FloatPerAtom, round_float_per_atom
from .multistage_opt import MultiStageOptMixin
from .workflow import MoleculeWorkflow

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

    @field_validator("constraints")
    def check_constraints(cls, constraints: Sequence[Constraint]) -> Sequence[Constraint]:
        if constraints:
            raise ValueError("ETKDG does not support constraints")

        return tuple(constraints)

    @field_validator("nci")
    def check_nci(cls, nci: bool) -> bool:
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


class ConformerGenMixin(BaseModel):
    """
    Mixin for workflows that need conformer generation.

    :param conf_gen_mode: Mode for calculations
    :param conf_gen_settings: settings for conformer generation
    :param constraints: constraints to add
    :param nci: add a constraining potential for non-covalent interactions
    :param max_confs: maximum number of conformers to keep
    """

    conf_gen_mode: Mode = Mode.RAPID
    conf_gen_settings: ConformerGenSettings = _sentinel  # type: ignore [assignment]
    constraints: Sequence[Constraint] = tuple()
    nci: bool = False
    max_confs: int | None = None

    @model_validator(mode="after")
    def validate_and_build_conf_gen_settings(self) -> Self:
        """Validate and build the ConformerGenSettings."""
        if self.conf_gen_settings is not _sentinel and self.conf_gen_mode != Mode.MANUAL:
            raise ValueError("Cannot specify conf_gen_settings with non-MANUAL mode")

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


class ConformerSearchWorkflow(ConformerSearchMixin, MoleculeWorkflow):
    """
    ConformerSearch Workflow.

    Inherited:
    :param initial_molecule: Molecule of interest
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

    # Results
    conformer_uuids: list[list[UUID | None]] = Field(default_factory=list)
    energies: Annotated[FloatPerAtom, AfterValidator(round_float_per_atom(6))] = Field(default_factory=list)
