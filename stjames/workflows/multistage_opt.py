from typing import Any, Sequence, TypeVar

from pydantic import model_validator

from ..constraint import Constraint
from ..method import XTB_METHODS, Method
from ..mode import Mode
from ..opt_settings import OptimizationSettings
from ..settings import Settings
from ..solvent import Solvent, SolventSettings
from ..task import Task
from .workflow import Workflow


class MultiStageOpt(Workflow):
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

    Note: allows a single point to be called with no optimization.

    Inherited
    :param initial_molecule: Molecule of interest

    :param mode: Mode for calculations
    :param optimization_settings: list of opt settings to apply successively
    :param singlepoint_settings: final single point settings
    :param solvent: solvent to use
    :param xtb_preopt: pre-optimize with xtb (sets based on mode when None)
    :param constraints: constraints for optimization
    :param transition_state: whether this is a transition state

    >>> from stjames.molecule import Molecule, Atom
    >>> He = Molecule(charge=0, multiplicity=1, atoms=[Atom(atomic_number=2, position=[0, 0, 0])])
    >>> mso = MultiStageOpt(initial_molecule=He, mode=Mode.RAPID, solvent="water")
    >>> mso.level_of_theory
    'r2scan_3c/cpcm(water)//gfn2_xtb'
    """

    mode: Mode | None = None
    optimization_settings: Sequence[Settings] = tuple()
    singlepoint_settings: Settings | None = None
    solvent: Solvent | None = None
    xtb_preopt: bool | None = None
    constraints: Sequence[Constraint] = tuple()
    transition_state: bool = False

    def __str__(self) -> str:
        if self.mode != Mode.MANUAL:
            return f"<MultiStageOpt {self.mode}>"

        return f"<MultiStageOpt {self.level_of_theory}>"

    @property
    def level_of_theory(self) -> str:
        methods = []
        if self.singlepoint_settings:
            methods.append(self.singlepoint_settings)
        methods += reversed(self.optimization_settings)

        return "//".join(m.level_of_theory for m in methods)

    @model_validator(mode="before")
    def check_only_mode_or_settings(cls, values: dict[str, Any]) -> dict[str, Any]:
        opt = values.get("optimization_settings")
        sp = values.get("singlepoint_settings")
        opt_or_sp_set = opt is not None or sp is not None

        match mode := values.get("mode"):
            case None:
                values["mode"] = Mode.MANUAL
                if not opt_or_sp_set:
                    raise ValueError("Must specify at least one of optimization_settings, singlepoint_settings, or mode")
            case Mode.MANUAL:
                if not opt_or_sp_set:
                    raise ValueError("Must specify at least one of optimization_settings or singlepoint_settings with MANUAL mode")
            case Mode.DEBUG:
                raise NotImplementedError(f"Unsupported mode: {mode}")
            case _:
                if opt_or_sp_set:
                    raise ValueError(f"Cannot specify optimization_settings or singlepoint_settings with {mode=}")
                values = _assign_settings_by_mode(values)

        return values


_T = TypeVar("_T", bound=dict[str, Any])


def _assign_settings_by_mode(values: _T) -> _T:
    """
    Construct the settings based on the mode.

    Defaults to RAPID mode.
    """
    assert values["mode"]
    if values["mode"] == Mode.MANUAL:
        return values

    opt_settings = OptimizationSettings()  # constraints=values.get("constraints"))

    # No solvent in pre-opt
    xtb_preopt = values.get("xtb_preopt")
    OPT = [Task.OPTIMIZE if not values.get("transition_state") else Task.OPTIMIZE_TS]
    gfn0_pre_opt = [Settings(method=Method.GFN0_XTB, tasks=OPT, mode=Mode.RAPID, opt_settings=opt_settings)]
    gfn2_pre_opt = [Settings(method=Method.GFN2_XTB, tasks=OPT, mode=Mode.RAPID, opt_settings=opt_settings)]

    def opt(method: Method, basis_set: str | None = None, solvent: Solvent | None = None, freq: bool = False) -> Settings:
        """Generates optimization settings."""
        model = "gbsa" if method in XTB_METHODS else "cpcm"

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
            opt_settings=opt_settings,
        )

    solvent = values.get("solvent")
    match mode := values["mode"]:
        case Mode.RECKLESS:
            # no-pre-opt
            optimization_settings = [opt(Method.GFN_FF)]
            singlepoint_settings = sp(Method.GFN2_XTB, solvent=solvent)

        case Mode.RAPID:
            optimization_settings = [
                *gfn0_pre_opt * bool(xtb_preopt),
                opt(Method.GFN2_XTB, freq=True),
            ]
            singlepoint_settings = sp(Method.R2SCAN3C, solvent=solvent)

        case Mode.CAREFUL:
            optimization_settings = [
                *gfn2_pre_opt * (xtb_preopt is None),
                opt(Method.B973C, solvent=solvent, freq=True),
            ]
            singlepoint_settings = sp(Method.WB97X3C, solvent=solvent)

        case Mode.METICULOUS:
            optimization_settings = [
                *gfn2_pre_opt * (xtb_preopt is None),
                opt(Method.B973C, solvent=solvent),
                opt(Method.WB97X3C, solvent=solvent, freq=True),
            ]
            singlepoint_settings = sp(Method.WB97MD3BJ, "def2-TZVPPD", solvent=solvent)

        case _:
            raise NotImplementedError(f"Unsupported {mode=}")

    values["optimization_settings"] = optimization_settings
    values["singlepoint_settings"] = singlepoint_settings

    return values
