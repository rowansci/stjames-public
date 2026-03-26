"""Engine/method/correction/task compatibility tables.

Matches the rules enforced in tinbergen's ``engines.tsx``. Used by
``Settings`` validators to catch invalid combinations at construction time.
"""

from .correction import Correction
from .engine import Engine
from .method import Method
from .task import Task

ENGINE_METHODS: dict[Engine, frozenset[Method]] = {
    Engine.AIMNET2: frozenset(
        {
            Method.AIMNET2_WB97MD3,
        }
    ),
    Engine.EGRET: frozenset(
        {
            Method.EGRET_1,
            Method.EGRET_1E,
            Method.EGRET_1T,
        }
    ),
    Engine.OMOL25: frozenset(
        {
            Method.OMOL25_CONSERVING_S,
            Method.UMA_S_OMOL,
            Method.UMA_S_1_2_OMOL,
            Method.UMA_M_OMOL,
            Method.UMA_S_OMAT,
            Method.UMA_S_1_2_OMAT,
            Method.UMA_M_OMAT,
            Method.UMA_S_OMC,
            Method.UMA_S_1_2_OMC,
            Method.UMA_M_OMC,
        }
    ),
    Engine.ORB: frozenset(
        {
            Method.ORB_V3_CONSERVATIVE_INF_OMAT,
            Method.ORB_V3_CONSERVATIVE_OMOL,
        }
    ),
    Engine.XTB: frozenset(
        {
            Method.GFN_FF,
            Method.GFN0_XTB,
            Method.GFN1_XTB,
            Method.GFN2_XTB,
            Method.G_XTB,
        }
    ),
    Engine.TBLITE: frozenset(
        {
            Method.GFN2_XTB,
        }
    ),
    Engine.OPENFF: frozenset(
        {
            Method.OFF_SAGE_2_0_0,
            Method.OFF_SAGE_2_2_1,
            Method.OFF_SAGE_2_3_0,
            Method.SMIRNOFF_2_0_0_AMBER_AM1BCC,
            Method.SMIRNOFF_2_2_1_AMBER_AM1BCC,
        }
    ),
    Engine.PSI4: frozenset(
        {
            Method.HARTREE_FOCK,
            Method.PBE,
            Method.R2SCAN,
            Method.TPSS,
            Method.M06L,
            Method.PBE0,
            Method.B3LYP,
            Method.TPSSH,
            Method.M06,
            Method.M062X,
            Method.CAMB3LYP,
            Method.WB97XV,
            Method.WB97XD3,
            Method.WB97MV,
            Method.WB97MD3BJ,
            Method.HF3C,
            Method.B973C,
            Method.B97D3BJ,
            Method.R2SCAN3C,
            Method.WB97X3C,
            Method.DSDBLYPD3BJ,
        }
    ),
    Engine.PYSCF: frozenset(
        {
            Method.HARTREE_FOCK,
            Method.PBE,
            Method.R2SCAN,
            Method.TPSS,
            Method.M06L,
            Method.PBE0,
            Method.B3LYP,
            Method.TPSSH,
            Method.M06,
            Method.M062X,
            Method.CAMB3LYP,
            Method.WB97XV,
            Method.WB97MV,
            Method.WB97MD3BJ,
            Method.SKALA,
        }
    ),
    Engine.GPU4PYSCF: frozenset(
        {
            Method.HARTREE_FOCK,
            Method.PBE,
            Method.R2SCAN,
            Method.TPSS,
            Method.M06L,
            Method.PBE0,
            Method.B3LYP,
            Method.TPSSH,
            Method.M06,
            Method.M062X,
            Method.CAMB3LYP,
            Method.WB97XV,
            Method.WB97MV,
            Method.WB97MD3BJ,
            Method.SKALA,
        }
    ),
}

# Engine-level defaults; use get_supported_corrections() for method-specific overrides.
ENGINE_CORRECTIONS: dict[Engine, frozenset[Correction]] = {
    Engine.AIMNET2: frozenset(),
    Engine.EGRET: frozenset(),
    Engine.OMOL25: frozenset(),
    Engine.ORB: frozenset({Correction.D3}),
    Engine.XTB: frozenset(),
    Engine.TBLITE: frozenset(),
    Engine.OPENFF: frozenset(),
    Engine.PSI4: frozenset({Correction.D3BJ}),
    Engine.PYSCF: frozenset({Correction.D3, Correction.D3BJ, Correction.D4}),
    Engine.GPU4PYSCF: frozenset({Correction.D3, Correction.D3BJ, Correction.D4}),
}

# Per-(method, engine) overrides on top of ENGINE_CORRECTIONS.
# r2scan+psi4 and b3lyp allow D4; M06 family on psi4 only allows D3.
_METHOD_ENGINE_CORRECTION_OVERRIDES: dict[tuple[Method, Engine], frozenset[Correction]] = {
    (Method.R2SCAN, Engine.PSI4): frozenset({Correction.D3BJ, Correction.D4}),
    (Method.B3LYP, Engine.PSI4): frozenset({Correction.D3BJ, Correction.D4}),
    (Method.B3LYP, Engine.PYSCF): frozenset({Correction.D3, Correction.D3BJ, Correction.D4}),
    (Method.B3LYP, Engine.GPU4PYSCF): frozenset({Correction.D3, Correction.D3BJ, Correction.D4}),
    (Method.M06, Engine.PSI4): frozenset({Correction.D3}),
    (Method.M06L, Engine.PSI4): frozenset({Correction.D3}),
    (Method.M062X, Engine.PSI4): frozenset({Correction.D3}),
}

# Tasks unavailable per engine (periodic mode may disable additional tasks).
ENGINE_DISABLED_TASKS: dict[Engine, frozenset[Task]] = {
    Engine.AIMNET2: frozenset({Task.DIPOLE, Task.SPIN_DENSITY}),
    Engine.EGRET: frozenset({Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY}),
    Engine.OMOL25: frozenset({Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY}),
    Engine.ORB: frozenset({Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY}),
    Engine.XTB: frozenset({Task.SPIN_DENSITY}),
    Engine.TBLITE: frozenset({Task.SPIN_DENSITY, Task.FREQUENCIES, Task.OPTIMIZE_TS}),
    Engine.OPENFF: frozenset(),
    Engine.PSI4: frozenset(),
    Engine.PYSCF: frozenset(),
    Engine.GPU4PYSCF: frozenset(),
}

# Engines that do not support open-shell (multiplicity > 1) calculations.
ENGINE_NO_OPEN_SHELL: frozenset[Engine] = frozenset(
    {
        Engine.AIMNET2,
        Engine.EGRET,
        Engine.ORB,
    }
)

# Engines that accept an explicit basis set.
ENGINE_SUPPORTS_BASIS_SET: frozenset[Engine] = frozenset(
    {
        Engine.PSI4,
        Engine.PYSCF,
        Engine.GPU4PYSCF,
    }
)


def get_supported_corrections(method: Method, engine: Engine) -> frozenset[Correction]:
    """Return corrections supported for a given method/engine combination.

    Applies per-method overrides from ``_METHOD_ENGINE_CORRECTION_OVERRIDES``
    on top of ``ENGINE_CORRECTIONS``.

    :param method: computational method
    :param engine: compute engine
    :returns: supported corrections, empty if none
    """
    if (method, engine) in _METHOD_ENGINE_CORRECTION_OVERRIDES:
        return _METHOD_ENGINE_CORRECTION_OVERRIDES[(method, engine)]
    return ENGINE_CORRECTIONS.get(engine, frozenset())
