from typing import Self, Sequence

from pydantic import BaseModel, Field, model_validator

from ..constraint import Constraint
from ..method import XTB_METHODS, Method
from ..mode import Mode
from ..opt_settings import OptimizationSettings
from ..settings import Settings
from ..solvent import Solvent, SolventSettings
from ..task import Task
from ..types import UUID
from .workflow import Workflow


class MultiStageOptSettings(BaseModel):
    """
    Settings for multi-stage optimizations.

    RECKLESS
        GFN2-xTB//GFN-FF (no pre-opt)
    RAPID *default
        r²SCAN-3c//GFN2-xTB with GFN0-xTB pre-opt (off by default)
    CAREFUL
        wB97X-3c//r²SCAN-3c with GFN2-xTB pre-opt
    METICULOUS
        wB97M-D3BJ/def2-TZVPPD//wB97X-3c//r²SCAN-3c with GFN2-xTB pre-opt

    Notes:
    - No solvent in any optimizations when using Modes
    - If solvent: xTB singlepoints use CPCMX, xTB optimizations use ALBP, all else use CPCM
    - Allows a single point to be called with no optimization

    :param mode: Mode for settings
    :param optimization_settings: list of opt settings to apply successively
    :param singlepoint_settings: final single point settings
    :param solvent: solvent to use for singlepoint
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param constraints: constraints for optimization
    :param transition_state: whether this is a transition state
    :param frequencies: whether to calculate frequencies

    >>> msos = MultiStageOptSettings(mode=Mode.RAPID, solvent="water")
    >>> msos
    <MultiStageOptSettings RAPID>
    >>> msos.level_of_theory
    'r2scan_3c/cpcm(water)//gfn2_xtb'
    """

    mode: Mode
    optimization_settings: Sequence[Settings] = tuple()
    singlepoint_settings: Settings | None = None
    solvent: Solvent | None = None
    xtb_preopt: bool = False
    constraints: Sequence[Constraint] = tuple()
    transition_state: bool = False
    frequencies: bool = False

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        """
        String representation of the settings.

        >>> print(MultiStageOptSettings(mode=Mode.RAPID, solvent="water"))
        <MultiStageOptSettings RAPID>
        """
        if self.mode != Mode.MANUAL:
            return f"<{type(self).__name__} {self.mode.name}>"

        return f"<{type(self).__name__} {self.level_of_theory}>"

    @property
    def level_of_theory(self) -> str:
        """
        Returns the level of theory for the workflow.

        >>> msos = MultiStageOptSettings(mode=Mode.RAPID, solvent="hexane")
        >>> msos.level_of_theory
        'r2scan_3c/cpcm(hexane)//gfn2_xtb'
        """
        methods = [self.singlepoint_settings] if self.singlepoint_settings else []
        methods += reversed(self.optimization_settings)

        return "//".join(m.level_of_theory for m in methods)

    @model_validator(mode="after")
    def set_mode_and_settings(self) -> Self:
        """Check mode and settings are properly specified, and assign settings based on mode."""
        opt_or_sp = bool(self.optimization_settings) or bool(self.singlepoint_settings)

        if self.mode == Mode.AUTO:
            self.mode = Mode.RAPID

        match (self.mode, opt_or_sp):
            case (Mode.DEBUG, _):
                raise NotImplementedError("Unsupported mode: DEBUG")

            case (Mode.MANUAL, False):
                raise ValueError("Must specify at least one of optimization_settings or singlepoint_settings with MANUAL mode")
            case (Mode.MANUAL, True):
                pass

            case (mode, True):
                raise ValueError(f"Cannot specify optimization_settings or singlepoint_settings with {mode=}")

            case (mode, False):
                self._assign_settings_by_mode(mode)

        return self

    def _assign_settings_by_mode(self, mode: Mode) -> None:
        """
        Construct the settings based on the Mode.

        :param mode: Mode to use
        """
        opt_settings = OptimizationSettings(constraints=self.constraints, transition_state=self.transition_state)

        # No solvent in pre-opt
        OPT = [Task.OPTIMIZE if not self.transition_state else Task.OPTIMIZE_TS]
        gfn0_pre_opt = [Settings(method=Method.GFN0_XTB, tasks=OPT, mode=Mode.RAPID, opt_settings=opt_settings)]
        gfn2_pre_opt = [Settings(method=Method.GFN2_XTB, tasks=OPT, mode=Mode.RAPID, opt_settings=opt_settings)]

        def opt(method: Method, basis_set: str | None = None, solvent: Solvent | None = None, freq: bool = False) -> Settings:
            """Generates optimization settings."""
            model = "alpb" if method in XTB_METHODS else "cpcm"

            return Settings(
                method=method,
                basis_set=basis_set,
                tasks=OPT + [Task.FREQUENCIES] * freq,
                solvent_settings=SolventSettings(solvent=solvent, model=model) if solvent else None,
                opt_settings=opt_settings,
            )

        def sp(method: Method, basis_set: str | None = None, solvent: Solvent | None = None) -> Settings:
            """Generate singlepoint settings."""
            model = "cpcmx" if method in XTB_METHODS else "cpcm"

            return Settings(
                method=method,
                basis_set=basis_set,
                solvent_settings=SolventSettings(solvent=solvent, model=model) if solvent else None,
            )

        match mode:
            case Mode.RECKLESS:
                self.xtb_preopt = False
                self.optimization_settings = [opt(Method.GFN_FF, freq=self.frequencies)]
                self.singlepoint_settings = sp(Method.GFN2_XTB, solvent=self.solvent)

            case Mode.RAPID:
                self.optimization_settings = [
                    *gfn0_pre_opt * self.xtb_preopt,
                    opt(Method.GFN2_XTB, freq=self.frequencies),
                ]
                self.singlepoint_settings = sp(Method.R2SCAN3C, solvent=self.solvent)

            case Mode.CAREFUL:
                self.optimization_settings = [
                    *gfn2_pre_opt * self.xtb_preopt,
                    opt(Method.R2SCAN3C, freq=self.frequencies),
                ]
                self.singlepoint_settings = sp(Method.WB97X3C, solvent=self.solvent)

            case Mode.METICULOUS:
                self.optimization_settings = [
                    *gfn2_pre_opt * self.xtb_preopt,
                    opt(Method.R2SCAN3C),
                    opt(Method.WB97X3C, freq=self.frequencies),
                ]
                self.singlepoint_settings = sp(Method.WB97MD3BJ, "def2-TZVPPD", solvent=self.solvent)

            case mode:
                raise NotImplementedError(f"Cannot assign settings for {mode=}")


class MultiStageOptWorkflow(Workflow, MultiStageOptSettings):
    """
    Workflow for multi-stage optimizations.

    Inherited
    :param initial_molecule: Molecule of interest
    :param mode: Mode for workflow
    :param optimization_settings: list of opt settings to apply successively
    :param singlepoint_settings: final single point settings
    :param solvent: solvent to use for singlepoint
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param constraints: constraints for optimization
    :param transition_state: whether this is a transition state
    :param frequencies: whether to calculate frequencies

    Populated while running
    :param calculations: list of calculation UUIDs

    >>> from stjames.molecule import Atom, Molecule
    >>> He = Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
    >>> msow = MultiStageOptWorkflow(initial_molecule=He, mode=Mode.RAPID, solvent="water")
    >>> msow
    <MultiStageOptWorkflow RAPID>
    >>> msow.level_of_theory
    'r2scan_3c/cpcm(water)//gfn2_xtb'
    """

    # Populated while running the workflow
    calculations: list[UUID | None] = Field(default_factory=list)

    def __repr__(self) -> str:
        if self.mode != Mode.MANUAL:
            return f"<{type(self).__name__} {self.mode.name}>"

        return f"<{type(self).__name__} {self.level_of_theory}>"


# the id of a mutable object may change, thus using object()
_sentinel_msos = object()


class MultiStageOptMixin(BaseModel):
    """
    Mixin for workflows that use MultiStageOptSettings.
    """

    mso_mode: Mode = Mode.AUTO
    # Need to use a sentinel object to make both mypy and pydantic happy
    multistage_opt_settings: MultiStageOptSettings = _sentinel_msos  # type: ignore [assignment]
    solvent: Solvent | None = None
    xtb_preopt: bool = False
    constraints: Sequence[Constraint] = tuple()
    transition_state: bool = False
    frequencies: bool = False

    @model_validator(mode="after")
    def set_mso_settings(self) -> Self:
        if self.mso_mode == Mode.AUTO:
            self.mso_mode = Mode.RAPID

        match self.mso_mode, self.multistage_opt_settings:
            case (Mode.DEBUG, _):
                raise NotImplementedError("Unsupported mode: DEBUG")

            case (Mode.MANUAL, msos) if msos is _sentinel_msos:
                raise ValueError("Must specify multistage_opt_settings with MANUAL mode")
            case (Mode.MANUAL, _):
                pass

            case (mso_mode, msos) if msos is not _sentinel_msos:
                raise ValueError(f"Cannot specify multistage_opt_settings with {mso_mode=}, {msos=}")

            case (mso_mode, _):
                self.multistage_opt_settings = MultiStageOptSettings(
                    mode=mso_mode,
                    solvent=self.solvent,
                    xtb_preopt=self.xtb_preopt,
                    constraints=self.constraints,
                    transition_state=self.transition_state,
                    frequencies=self.frequencies,
                )

        return self


def build_mso_settings(
    sp_method: Method,
    sp_basis_set: str | None,
    opt_methods: list[Method],
    opt_basis_sets: list[str | None],
    mode: Mode = Mode.MANUAL,
    solvent: Solvent | None = None,
    use_solvent_for_opt: bool = False,
    constraints: list[Constraint] | None = None,
    transition_state: bool = False,
    frequencies: bool = False,
) -> MultiStageOptSettings:
    """
    Helper function to construct multi-stage opt settings objects manually.

    There's no xTB pre-optimization here - add that yourself!

    :param optimization_settings: optimization settings to apply successively
    :param singlepoint_settings: final single point settings
    :param mode: Mode for settings, defaults to `MANUAL`
    :param solvent: solvent to use
    :param use_solvent_for_opt: whether to conduct opts with solvent
    :param constraints: constraints for optimization
    :param transition_state: whether this is a transition state
    :param frequencies: whether to calculate frequencies
    :returns: MultiStageOptSettings
    """
    if constraints is None:
        constraints = []

    opt_settings = OptimizationSettings(constraints=constraints, transition_state=transition_state)

    OPT = [Task.OPTIMIZE if not transition_state else Task.OPTIMIZE_TS]

    def opt(method: Method, basis_set: str | None = None, solvent: Solvent | None = None, freq: bool = False) -> Settings:
        """Generates optimization settings."""
        model = "alpb" if method in XTB_METHODS else "cpcm"

        return Settings(
            method=method,
            basis_set=basis_set,
            tasks=OPT + [Task.FREQUENCIES] * freq,
            solvent_settings=SolventSettings(solvent=solvent, model=model) if (solvent and use_solvent_for_opt) else None,
            opt_settings=opt_settings,
        )

    def sp(method: Method, basis_set: str | None = None, solvent: Solvent | None = None) -> Settings:
        """Generate singlepoint settings."""
        model = "cpcmx" if method in XTB_METHODS else "cpcm"

        return Settings(
            method=method,
            basis_set=basis_set,
            tasks=[Task.ENERGY],
            solvent_settings=SolventSettings(solvent=solvent, model=model) if solvent else None,
        )

    return MultiStageOptSettings(
        mode=mode,
        optimization_settings=[
            opt(method=method, basis_set=basis_set, solvent=solvent, freq=frequencies) for method, basis_set in zip(opt_methods, opt_basis_sets, strict=True)
        ],
        singlepoint_settings=sp(method=sp_method, basis_set=sp_basis_set, solvent=solvent),
        solvent=solvent,
        xtb_preopt=False,
        constraints=constraints,
        transition_state=transition_state,
        frequencies=frequencies,
    )


def multi_stage_opt_settings_from_workflow(msow: MultiStageOptWorkflow) -> MultiStageOptSettings:
    """
    Helper function to convert a MultiStageOptWorkflow to MultiStageOptSettings.

    :param msow: MultiStageOptWorkflow
    :returns: MultiStageOptSettings

    >>> from stjames.molecule import Atom, Molecule
    >>> He = Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
    >>> msos = multi_stage_opt_settings_from_workflow(
    ...     MultiStageOptWorkflow(initial_molecule=He, mode=Mode.RAPID, solvent="water")
    ... )
    >>> print(msos)
    <MultiStageOptSettings RAPID>
    >>> msos.level_of_theory
    'r2scan_3c/cpcm(water)//gfn2_xtb'
    >>> msos.solvent
    <Solvent.WATER: 'water'>
    >>> msos.xtb_preopt
    False
    """
    data = dict(msow)
    del data["calculations"]
    return MultiStageOptSettings.construct(**data)
