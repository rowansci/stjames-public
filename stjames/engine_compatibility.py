"""Engine/method/correction/task compatibility tables.

Used by ``Settings`` validators to catch invalid combinations at construction time.
"""

from itertools import product

from .correction import Correction
from .engine import Engine
from .method import Method
from .solvent import SolventModel
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
    Engine.GPU4PYSCF: frozenset(
        {
            Method.HARTREE_FOCK,
            Method.BP86,
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
            Method.WB97X3C,
            Method.SKALA,
        }
    ),
    Engine.MOPAC: frozenset(
        {
            Method.PM6,
            Method.PM6_D3H4X,
            Method.PM6_ORG,
            Method.PM7,
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
            Method.BP86,
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
            Method.WB97X3C,
            Method.SKALA,
        }
    ),
    Engine.QUANTUM_ESPRESSO: frozenset(
        {
            Method.HARTREE_FOCK,
            Method.PBE,
            Method.BP86,
            Method.R2SCAN,
            Method.TPSS,
            Method.M06L,
            Method.B97D3BJ,
            Method.PBE0,
            Method.B3LYP,
            Method.TPSSH,
            Method.M06,
            Method.M062X,
            Method.CAMB3LYP,
            Method.WB97XD3,
        }
    ),
    Engine.TBLITE: frozenset(
        {
            Method.GFN1_XTB,
            Method.GFN2_XTB,
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
}

METHOD_ENGINES: dict[Method, frozenset[Engine]] = {
    method: frozenset(engine for engine, methods in ENGINE_METHODS.items() if method in methods)
    for method in {m for methods in ENGINE_METHODS.values() for m in methods}
}

# Engine-level defaults; use get_supported_corrections() for method-specific overrides.
ENGINE_CORRECTIONS: dict[Engine, frozenset[Correction]] = {
    Engine.AIMNET2: frozenset(),
    Engine.EGRET: frozenset(),
    Engine.GPU4PYSCF: frozenset({Correction.D3, Correction.D3BJ, Correction.D4}),
    Engine.MOPAC: frozenset(),
    Engine.OMOL25: frozenset(),
    Engine.OPENFF: frozenset(),
    Engine.ORB: frozenset({Correction.D3}),
    Engine.PSI4: frozenset({Correction.D3BJ}),
    Engine.PYSCF: frozenset({Correction.D3, Correction.D3BJ, Correction.D4}),
    Engine.QUANTUM_ESPRESSO: frozenset({Correction.D3, Correction.D3BJ}),
    Engine.TBLITE: frozenset(),
    Engine.XTB: frozenset(),
}

# Per-(method, engine) overrides on top of ENGINE_CORRECTIONS.
# r2scan+psi4 and b3lyp allow D4; M06 family on psi4 only allows D3.
# Methods with dispersion baked in (WB97X-D3, WB97X-3C, WB97M-D3BJ, etc.) disallow further corrections.
_METHOD_ENGINE_CORRECTION_OVERRIDES: dict[tuple[Method, Engine], frozenset[Correction]] = {
    (Method.R2SCAN, Engine.PSI4): frozenset({Correction.D3BJ, Correction.D4}),
    (Method.B3LYP, Engine.PSI4): frozenset({Correction.D3BJ, Correction.D4}),
    (Method.B3LYP, Engine.PYSCF): frozenset({Correction.D3, Correction.D3BJ, Correction.D4}),
    (Method.B3LYP, Engine.GPU4PYSCF): frozenset({Correction.D3, Correction.D3BJ, Correction.D4}),
    (Method.M06, Engine.PSI4): frozenset({Correction.D3}),
    (Method.M06L, Engine.PSI4): frozenset({Correction.D3}),
    (Method.M062X, Engine.PSI4): frozenset({Correction.D3}),
    # Dispersion baked in — no additional corrections allowed
    **{
        (method, engine): frozenset()
        for method, engine in product(
            (Method.WB97XD3, Method.WB97X3C, Method.WB97MD3BJ, Method.WB97XV, Method.WB97MV, Method.DSDBLYPD3BJ, Method.B97D3BJ),
            (Engine.PSI4, Engine.PYSCF, Engine.GPU4PYSCF),
        )
    },
}

# Tasks unavailable per engine (periodic mode may disable additional tasks).
ENGINE_DISABLED_TASKS: dict[Engine, frozenset[Task]] = {
    Engine.AIMNET2: frozenset({Task.DIPOLE, Task.SPIN_DENSITY}),
    Engine.EGRET: frozenset({Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY}),
    Engine.GPU4PYSCF: frozenset(),
    Engine.MOPAC: frozenset({Task.DIPOLE, Task.SPIN_DENSITY}),
    Engine.OMOL25: frozenset({Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY}),
    Engine.ORB: frozenset({Task.CHARGE, Task.DIPOLE, Task.SPIN_DENSITY}),
    Engine.OPENFF: frozenset(),
    Engine.PSI4: frozenset(),
    Engine.PYSCF: frozenset(),
    Engine.QUANTUM_ESPRESSO: frozenset({Task.DIPOLE, Task.SPIN_DENSITY, Task.HESSIAN, Task.FREQUENCIES}),
    Engine.TBLITE: frozenset({Task.SPIN_DENSITY}),
    Engine.XTB: frozenset({Task.SPIN_DENSITY}),
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
        Engine.GPU4PYSCF,
        Engine.PSI4,
        Engine.PYSCF,
        Engine.QUANTUM_ESPRESSO,
    }
)


# Solvent models supported per engine. Engines not listed do not support solvent models.
ENGINE_SOLVENT_MODELS: dict[Engine, frozenset[SolventModel]] = {
    Engine.AIMNET2: frozenset({SolventModel.ALPB, SolventModel.CPCMX}),
    Engine.EGRET: frozenset({SolventModel.ALPB, SolventModel.CPCMX}),
    Engine.GPU4PYSCF: frozenset({SolventModel.CPCM, SolventModel.PCM}),
    Engine.MOPAC: frozenset({SolventModel.COSMO}),
    Engine.OMOL25: frozenset({SolventModel.ALPB, SolventModel.CPCMX}),
    Engine.ORB: frozenset({SolventModel.ALPB, SolventModel.CPCMX}),
    Engine.PSI4: frozenset({SolventModel.COSMO, SolventModel.CPCM, SolventModel.PCM}),
    Engine.PYSCF: frozenset({SolventModel.COSMO, SolventModel.CPCM, SolventModel.PCM}),
    Engine.XTB: frozenset({SolventModel.ALPB, SolventModel.GBSA, SolventModel.CPCMX}),
}


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
