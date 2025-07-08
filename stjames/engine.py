from .base import LowercaseStrEnum


class Mode(LowercaseStrEnum):
    AIMNET2 = "aimnet2"
    MACE = "mace"
    OCP24 = "ocp24"
    OMOL25 = "omol25"
    ORB = "orb"
    TBLITE = "tblite"
    XTB = "xtb"
    TERRACHEM = "terrachem"
    PYSCF = "pyscf"
    PSI4 = "psi4"
    OPENFF = "openff"
    EGRET = "egret"
