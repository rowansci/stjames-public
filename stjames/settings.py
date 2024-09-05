from typing import Any, Optional, TypeVar

import pydantic

from .base import Base, UniqueList
from .basis_set import BasisSet
from .correction import Correction
from .method import METHODS_WITH_CORRECTION, PREPACKAGED_METHODS, Method
from .mode import Mode
from .opt_settings import OptimizationSettings
from .scf_settings import SCFSettings
from .solvent import SolventSettings
from .task import Task
from .thermochem_settings import ThermochemistrySettings

_T = TypeVar("_T")


class Settings(Base):
    method: Method = Method.HARTREE_FOCK
    basis_set: Optional[BasisSet] = None
    tasks: UniqueList[Task] = [Task.ENERGY, Task.CHARGE, Task.DIPOLE]
    corrections: UniqueList[Correction] = []

    mode: Mode = Mode.AUTO

    solvent_settings: Optional[SolventSettings] = None

    # scf/opt settings will be set automatically based on mode, but can be overridden manually
    scf_settings: SCFSettings = SCFSettings()
    opt_settings: OptimizationSettings = OptimizationSettings()
    thermochem_settings: ThermochemistrySettings = ThermochemistrySettings()

    # mypy has this dead wrong (https://docs.pydantic.dev/2.0/usage/computed_fields/)
    # Python 3.12 narrows the reason for the ignore to prop-decorator
    @pydantic.computed_field  # type: ignore[misc, prop-decorator, unused-ignore]
    @property
    def level_of_theory(self) -> str:
        corrections = list(filter(lambda x: x not in (None, ""), self.corrections))

        if self.method in PREPACKAGED_METHODS or self.basis_set is None:
            method = self.method.value
        elif self.method in METHODS_WITH_CORRECTION or len(corrections) == 0:
            method = f"{self.method.value}/{self.basis_set.name.lower()}"
        else:
            method = f"{self.method.value}-{'-'.join([c.value for c in corrections])}/{self.basis_set.name.lower()}"

        if self.solvent_settings is not None:
            method += f"/{self.solvent_settings.model.value}({self.solvent_settings.solvent.value})"

        return method

    def model_post_init(self, __context: Any) -> None:
        _assign_settings_by_mode(self)

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

    @pydantic.field_validator("basis_set", mode="before")
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

    @pydantic.field_validator("corrections", mode="before")
    @classmethod
    def remove_empty_string(cls, v: list[_T]) -> list[_T]:
        """Remove empty string values."""
        return [c for c in v if c] if v is not None else v


def _assign_settings_by_mode(settings: Settings) -> None:
    """Modifies ``scf_settings`` and ``opt_settings`` based on preset ``mode``."""
    mode = settings.mode

    if mode == Mode.AUTO:
        if (Task.OPTIMIZE in settings.tasks) or (Task.GRADIENT in settings.tasks) or (Task.FREQUENCIES in settings.tasks) or (Task.HESSIAN in settings.tasks):
            # noisy gradient! struggles to converge
            if settings.method == Method.AIMNET2_WB97MD3:
                mode = Mode.RAPID
            else:
                mode = Mode.CAREFUL
        else:
            mode = Mode.RAPID
    elif mode == Mode.MANUAL:
        return

    # modify scf settings!
    #
    # values based off of the following sources:
    # qchem:
    #   https://manual.q-chem.com/5.2/Ch4.S3.SS2.html
    #   https://manual.q-chem.com/5.2/Ch4.S5.SS2.html
    #
    # gaussian:
    #   https://gaussian.com/integral/
    #   https://gaussian.com/overlay5/
    #
    # orca:
    #   manual 4.2.1, §9.6.1 and §9.7.3
    #
    # psi4:
    #   https://psicode.org/psi4manual/master/autodir_options_c/module__scf.html
    #   https://psicode.org/psi4manual/master/autodoc_glossary_options_c.html
    #
    # terachem:
    #   manual, it's easy to locate everything.
    #
    # the below values are my best attempt at homogenizing various sources.
    # in general, eri_threshold should be 3 OOM lower than scf convergence
    scf_settings = settings.scf_settings
    if mode == Mode.RECKLESS:
        scf_settings.energy_threshold = 1e-5
        scf_settings.rms_error_threshold = 1e-7
        scf_settings.max_error_threshold = 1e-5
        scf_settings.rebuild_frequency = 100
        scf_settings.int_settings.eri_threshold = 1e-8
        scf_settings.int_settings.csam_multiplier = 3.0
        scf_settings.int_settings.pair_overlap_threshold = 1e-8
    elif mode == Mode.RAPID:
        scf_settings.energy_threshold = 5e-5
        scf_settings.rms_error_threshold = 1e-8
        scf_settings.max_error_threshold = 1e-6
        scf_settings.rebuild_frequency = 20
        scf_settings.int_settings.eri_threshold = 1e-9
        scf_settings.int_settings.csam_multiplier = 1.0
        scf_settings.int_settings.pair_overlap_threshold = 1e-9
    elif mode == Mode.CAREFUL:
        scf_settings.energy_threshold = 1e-6
        scf_settings.rms_error_threshold = 1e-9
        scf_settings.max_error_threshold = 1e-7
        scf_settings.rebuild_frequency = 10
        scf_settings.int_settings.eri_threshold = 1e-10
        scf_settings.int_settings.csam_multiplier = 1.0
        scf_settings.int_settings.pair_overlap_threshold = 1e-10
    elif mode == Mode.METICULOUS:
        scf_settings.energy_threshold = 1e-8
        scf_settings.rms_error_threshold = 1e-9
        scf_settings.max_error_threshold = 1e-7
        scf_settings.rebuild_frequency = 5
        scf_settings.int_settings.eri_threshold = 1e-12
        scf_settings.int_settings.csam_multiplier = 1.0
        scf_settings.int_settings.pair_overlap_threshold = 1e-12
    elif mode == Mode.DEBUG:
        scf_settings.energy_threshold = 1e-9
        scf_settings.rms_error_threshold = 1e-10
        scf_settings.max_error_threshold = 1e-9
        scf_settings.rebuild_frequency = 1
        scf_settings.int_settings.eri_threshold = 1e-14
        scf_settings.int_settings.csam_multiplier = 1e10  # in other words, disable CSAM
        scf_settings.int_settings.pair_overlap_threshold = 1e-14
    else:
        raise ValueError(f"Unknown mode ``{mode.value}``!")

    opt_settings = settings.opt_settings

    # constrained optimizations warrant loosening the settings a bit
    has_constraints = len(opt_settings.constraints) > 0

    # cf. DLFIND manual, and https://www.cup.uni-muenchen.de/ch/compchem/geom/basic.html
    # and the discussion at https://geometric.readthedocs.io/en/latest/how-it-works.html
    # in periodic systems, "normal" is 0.05 eV/Å ~= 2e-3 Hartree/Å, and "careful" is 0.01 ~= 4e-4
    if mode == Mode.RECKLESS:
        opt_settings.energy_threshold = 2e-5
        opt_settings.max_gradient_threshold = 7e-3
        opt_settings.rms_gradient_threshold = 6e-3
    elif mode == Mode.RAPID or (mode == Mode.CAREFUL and has_constraints):
        opt_settings.energy_threshold = 5e-5
        opt_settings.max_gradient_threshold = 5e-3
        opt_settings.rms_gradient_threshold = 3.5e-3
    elif mode == Mode.CAREFUL or (mode == Mode.METICULOUS and has_constraints):
        opt_settings.energy_threshold = 1e-6
        opt_settings.max_gradient_threshold = 9e-4
        opt_settings.rms_gradient_threshold = 6e-4
    elif mode == Mode.METICULOUS:
        opt_settings.energy_threshold = 1e-6
        opt_settings.max_gradient_threshold = 3e-5
        opt_settings.rms_gradient_threshold = 2e-5
    elif mode == Mode.DEBUG:
        opt_settings.energy_threshold = 1e-6
        opt_settings.max_gradient_threshold = 4e-6
        opt_settings.rms_gradient_threshold = 2e-6
    else:
        raise ValueError(f"Unknown mode ``{mode.value}``!")
