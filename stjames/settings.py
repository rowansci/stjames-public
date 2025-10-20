from typing import Any, Optional, Self, TypeVar

from pydantic import computed_field, field_validator, model_validator

from .base import Base, UniqueList
from .basis_set import BasisSet
from .compute_settings import ComputeSettings
from .correction import Correction
from .method import CORRECTABLE_NNP_METHODS, METHODS_WITH_CORRECTION, PREPACKAGED_METHODS, Method
from .mode import Mode
from .opt_settings import OptimizationSettings
from .scf_settings import SCFSettings
from .solvent import SolventSettings
from .task import Task
from .thermochem_settings import ThermochemistrySettings

_T = TypeVar("_T")


class Settings(Base):
    mode: Mode = Mode.AUTO

    method: Method = Method.HARTREE_FOCK
    basis_set: Optional[BasisSet] = None
    tasks: UniqueList[Task] = [Task.ENERGY, Task.CHARGE, Task.DIPOLE]
    corrections: UniqueList[Correction] = []

    solvent_settings: Optional[SolventSettings] = None

    # scf/opt settings will be set automatically based on mode, but can be overridden manually
    scf_settings: SCFSettings = SCFSettings()
    opt_settings: OptimizationSettings = OptimizationSettings()
    thermochem_settings: ThermochemistrySettings = ThermochemistrySettings()
    compute_settings: ComputeSettings = ComputeSettings()

    # mypy has this dead wrong (https://docs.pydantic.dev/2.0/usage/computed_fields/)
    # Python 3.12 narrows the reason for the ignore to prop-decorator
    @computed_field  # type: ignore[misc, prop-decorator, unused-ignore]
    @property
    def level_of_theory(self) -> str:
        corrections = list(filter(lambda x: x not in (None, ""), self.corrections))

        if self.method in CORRECTABLE_NNP_METHODS:
            method = self.method.value if not corrections else f"{self.method.value}-{'-'.join(c.value for c in corrections)}"
        elif self.method in PREPACKAGED_METHODS or self.basis_set is None:
            method = self.method.value
        elif self.method in METHODS_WITH_CORRECTION or not corrections:
            method = f"{self.method.value}/{self.basis_set.name.lower()}"
        else:
            method = f"{self.method.value}-{'-'.join(c.value for c in corrections)}/{self.basis_set.name.lower()}"

        if self.solvent_settings is not None:
            method += f"/{self.solvent_settings.model.value}({self.solvent_settings.solvent.value})"

        return method

    @field_validator("mode")
    @classmethod
    def set_mode_auto(cls, mode: Mode) -> Mode:
        """Set the mode to RAPID if AUTO is selected."""
        if mode == Mode.AUTO:
            return Mode.RAPID

        return mode

    @model_validator(mode="after")
    def validate_and_build(self) -> Self:
        if self.mode == Mode.AUTO:
            self.mode = Mode.RAPID

        self.opt_settings = _assign_opt_settings_by_mode(self.mode, self.opt_settings)

        return self

    def model_post_init(self, __context: Any) -> None:
        # figure out `optimize_ts`
        if Task.OPTIMIZE_TS in self.tasks:
            self.tasks.pop(self.tasks.index(Task.OPTIMIZE_TS))
            self.tasks.append(Task.OPTIMIZE)
            self.opt_settings.transition_state = True

        # composite methods have their own basis sets, so overwrite user stuff
        if self.method == Method.HF3C:
            self.basis_set = BasisSet(name="minix")
        elif self.method == Method.B973C:
            self.basis_set = BasisSet(name="def2-mTZVP")
        elif self.method == Method.R2SCAN3C:
            self.basis_set = BasisSet(name="def2-mTZVPP")
        elif self.method == Method.WB97X3C:
            self.basis_set = BasisSet(name="vDZP")

    @field_validator("basis_set", mode="before")
    @classmethod
    def parse_basis_set(cls, v: Any) -> BasisSet | dict[str, Any] | None:
        """Turn a string into a ``BasisSet`` object. (This is a little crude.)"""
        if isinstance(v, BasisSet):
            return None if v.name is None else v
        elif isinstance(v, dict):
            return None if v.get("name") is None else v
        elif isinstance(v, str):
            if len(v):
                return BasisSet(name=v)
            # "" is basically None, let's be real here...
            return None
        elif v is None:
            return None
        else:
            raise ValueError(f"invalid value ``{v}`` for ``basis_set``")

    @field_validator("corrections", mode="before")
    @classmethod
    def remove_empty_string(cls, v: list[_T]) -> list[_T]:
        """Remove empty string values."""
        return [c for c in v if c] if v is not None else v


def _assign_opt_settings_by_mode(mode: Mode, opt_settings: OptimizationSettings) -> OptimizationSettings:
    """
    Assign optimization settings based on the mode.

    Constraints lead to a lot of noise, so we need to loosen the thresholds.

    cf. DLFIND manual, and https://www.cup.uni-muenchen.de/ch/compchem/geom/basic.html
    and the discussion at https://geometric.readthedocs.io/en/latest/how-it-works.html
    in periodic systems, "normal" is 0.05 eV/Å ~= 2e-3 Hartree/Å, and "careful" is 0.01 ~= 4e-4

    Note: thresholds here are in units of Hartree/Å, not Hartree/Bohr as listed in many places.
    """
    opt_settings.energy_threshold = 1e-6
    match mode:
        case Mode.RECKLESS:
            opt_settings.energy_threshold = 2e-5
            opt_settings.max_gradient_threshold = 7e-3
            opt_settings.rms_gradient_threshold = 6e-3
        case Mode.RAPID:
            opt_settings.energy_threshold = 5e-5
            opt_settings.max_gradient_threshold = 5e-3
            opt_settings.rms_gradient_threshold = 3.5e-3
        case Mode.CAREFUL:
            opt_settings.max_gradient_threshold = 9e-4
            opt_settings.rms_gradient_threshold = 6e-4
        case Mode.METICULOUS:
            opt_settings.max_gradient_threshold = 3e-5
            opt_settings.rms_gradient_threshold = 2e-5
        case Mode.DEBUG:
            opt_settings.max_gradient_threshold = 4e-6
            opt_settings.rms_gradient_threshold = 2e-6
        case _:
            raise ValueError(f"Unknown mode ``{mode.value}``!")

    return opt_settings
