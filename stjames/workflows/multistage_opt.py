from typing import Self, Sequence

from pydantic import model_validator

from ..constraint import Constraint
from ..method import XTB_METHODS, Method
from ..mode import Mode
from ..opt_settings import OptimizationSettings
from ..settings import Settings
from ..solvent import Solvent, SolventSettings
from ..task import Task
from .workflow import UUID, Workflow


class MultiStageOptWorkflow(Workflow):
    """
    Workflow for multi-stage optimizations.

    RECKLESS
        GFN2-xTB//GFN-FF (no pre-opt)
    RAPID *default
        r²SCAN-3c//GFN2-xTB with GFN0-xTB pre-opt (off by default)
    CAREFUL
        wB97X-3c//B97-3c with GFN2-xTB pre-opt
    METICULOUS
        wB97M-D3BJ/def2-TZVPPD//wB97X-3c//B97-3c with GFN2-xTB pre-opt

    Notes:
    - No solvent in pre-opt
    - For solvent: xTB singlepoints use CPCMX, xTB optimizations use ALBP, all else use CPCM
    - Allows a single point to be called with no optimization

    Inherited
    :param initial_molecule: Molecule of interest

    :param mode: Mode for workflow
    :param optimization_settings: list of opt settings to apply successively
    :param singlepoint_settings: final single point settings
    :param solvent: solvent to use
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param constraints: constraints for optimization
    :param transition_state: whether this is a transition state

    >>> from stjames.molecule import Atom, Molecule
    >>> He = Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
    >>> mso = MultiStageOptWorkflow(initial_molecule=He, mode=Mode.RAPID, solvent="water")
    >>> mso
    <MultiStageOptWorkflow RAPID>
    >>> mso.level_of_theory
    'r2scan_3c/cpcm(water)//gfn2_xtb/alpb(water)'
    """

    mode: Mode
    optimization_settings: Sequence[Settings] = tuple()
    singlepoint_settings: Settings | None = None
    solvent: Solvent | None = None
    xtb_preopt: bool | None = None
    constraints: Sequence[Constraint] = tuple()
    transition_state: bool = False

    # Populated while running the workflow
    calculations: list[UUID] | None = None

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        """
        String representation of the workflow.

        >>> from stjames.molecule import Atom, Molecule
        >>> He = Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
        >>> print(MultiStageOptWorkflow(initial_molecule=He, mode=Mode.RAPID, solvent="water"))
        <MultiStageOptWorkflow RAPID>
        """
        if self.mode != Mode.MANUAL:
            return f"<MultiStageOptWorkflow {self.mode.name}>"

        return f"<MultiStageOptWorkflow {self.level_of_theory}>"

    @property
    def level_of_theory(self) -> str:
        """
        Returns the level of theory for the workflow.

        >>> from stjames.molecule import Atom, Molecule
        >>> He = Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
        >>> mso = MultiStageOptWorkflow(initial_molecule=He, mode=Mode.RAPID, solvent="hexane")
        >>> mso.level_of_theory
        'r2scan_3c/cpcm(hexane)//gfn2_xtb/alpb(hexane)'
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
                # no-pre-opt
                self.xtb_preopt = False
                self.optimization_settings = [opt(Method.GFN_FF, solvent=self.solvent)]
                self.singlepoint_settings = sp(Method.GFN2_XTB, solvent=self.solvent)

            case Mode.RAPID:
                self.xtb_preopt = bool(self.xtb_preopt)
                self.optimization_settings = [
                    *gfn0_pre_opt * self.xtb_preopt,
                    opt(Method.GFN2_XTB, freq=True, solvent=self.solvent),
                ]
                self.singlepoint_settings = sp(Method.R2SCAN3C, solvent=self.solvent)

            case Mode.CAREFUL:
                self.xtb_preopt = (self.xtb_preopt is None) or self.xtb_preopt
                self.optimization_settings = [
                    *gfn2_pre_opt * self.xtb_preopt,
                    opt(Method.B973C, solvent=self.solvent, freq=True),
                ]
                self.singlepoint_settings = sp(Method.WB97X3C, solvent=self.solvent)

            case Mode.METICULOUS:
                self.xtb_preopt = (self.xtb_preopt is None) or self.xtb_preopt
                self.optimization_settings = [
                    *gfn2_pre_opt * self.xtb_preopt,
                    opt(Method.B973C, solvent=self.solvent),
                    opt(Method.WB97X3C, solvent=self.solvent, freq=True),
                ]
                self.singlepoint_settings = sp(Method.WB97MD3BJ, "def2-TZVPPD", solvent=self.solvent)

            case mode:
                raise NotImplementedError(f"Cannot assign settings for {mode=}")

        assert self.xtb_preopt is not None
