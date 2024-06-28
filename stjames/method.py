from .base import LowercaseStrEnum


class Method(LowercaseStrEnum):
    HARTREE_FOCK = "hf"
    HF3C = "hf_3c"

    PBE = "pbe"
    B973C = "b97_3c"
    R2SCAN = "r2scan"
    R2SCAN3C = "r2scan3c"
    TPSS = "tpss"
    M06L = "m06l"

    PBE0 = "pbe0"
    B3LYP = "b3lyp"
    TPSSH = "tpssh"
    M06 = "m06"
    M062X = "m062x"

    CAMB3LYP = "camb3lyp"
    WB97XD3 = "wb97x_d3"
    WB97XV = "wb97x_v"
    WB97MV = "wb97m_v"
    WB97MD3BJ = "wb97m_d3bj"
    WB97X3C = "wb97x_3c"

    DSDBLYPD3BJ = "dsd_b3lyp_d3bj"

    AIMNET2_WB97MD3 = "aimnet2_wb97md3"

    # GFN0_XTB = "gfn1_xtb"
    GFN1_XTB = "gfn1_xtb"
    GFN2_XTB = "gfn2_xtb"
    GFN_FF = "gfn_ff"
