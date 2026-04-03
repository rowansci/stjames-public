from .base import LowercaseStrEnum


class Engine(LowercaseStrEnum):
    """Computational chemistry engine."""

    AIMNET2 = "aimnet2"
    EGRET = "egret"
    GPU4PYSCF = "gpu4pyscf"
    MACE = "mace"  # Deprecated
    MOPAC = "mopac"
    OPENFF = "openff"
    OMOL25 = "omol25"
    ORB = "orb"
    PSI4 = "psi4"
    PYSCF = "pyscf"
    QUANTUM_ESPRESSO = "quantum_espresso"
    TBLITE = "tblite"
    TERACHEM = "terachem"
    XTB = "xtb"
