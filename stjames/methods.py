from .base import LowercaseStrEnum


class Method(LowercaseStrEnum):
    #### Hartree–Fock:
    HARTREE_FOCK = "hf"
    HF3C = "hf-3c"

    #### DFT:
    # LDA
    LSDA = "lsda"

    # local GGA
    PBE = "pbe"
    BLYP = "blyp"
    BP86 = "bp86"
    B97D3 = "b97-d3"
    B973C = "b97-3c"

    # mGGA
    R2SCAN = "r2scan"
    TPSS = "tpss"
    M06L = "m06l"

    # hybrid GGA
    PBE0 = "pbe0"
    B3LYP = "b3lyp"
    B3PW91 = "b3pw91"
    B97MV = "b97m_v"

    # hybrid mGGA
    TPSS0 = "tpss0"
    M06 = "m06"
    M062X = "m062x"

    # range-separated
    CAMB3LYP = "camb3lyp"
    WB97XD = "wb97x_d"
    WB97XD3 = "wb97x_d3"
    WB97XV = "wb97x_v"
    WB97MV = "wb97m_v"

    # ML methods
    AIMNET2_WB97MD3 = "aimnet2_wb97md3"
    AIMNET2_B973C = "aimnet2_b973c"

    # xTB methods
    GFN1_XTB = "gfn1_xtb"
    GFN2_XTB = "gfn2_xtb"
