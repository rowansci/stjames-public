from .base import LowercaseStrEnum


class Engine(LowercaseStrEnum):
    AIMNET2 = "aimnet2"
    MACE = "mace"
    OMOL25 = "omol25"
    ORB = "orb"
    TBLITE = "tblite"
    XTB = "xtb"
    TERACHEM = "terachem"
    PYSCF = "pyscf"
    PSI4 = "psi4"
    OPENFF = "openff"
    EGRET = "egret"
