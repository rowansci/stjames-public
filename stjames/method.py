from typing import Literal

from .base import LowercaseStrEnum


class Method(LowercaseStrEnum):
    HARTREE_FOCK = "hf"
    HF3C = "hf_3c"

    PBE = "pbe"
    B973C = "b97_3c"
    B97D3BJ = "b97_d3bj"
    R2SCAN = "r2scan"
    R2SCAN3C = "r2scan_3c"
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

    DSDBLYPD3BJ = "dsd_blyp_d3bj"

    AIMNET2_WB97MD3 = "aimnet2_wb97md3"
    MACE_MP_0 = "mace_mp_0"
    OCP24_S = "ocp24_s"
    OCP24_L = "ocp24_l"

    GFN_FF = "gfn_ff"
    GFN0_XTB = "gfn0_xtb"
    GFN1_XTB = "gfn1_xtb"
    GFN2_XTB = "gfn2_xtb"

    # this was going to be removed, but Jonathon wrote such a nice basis set test... it's off the front end.
    BP86 = "bp86"

    OFF_SAGE_2_2_1 = "off_sage_2_2_1"


NNPMethod = Literal[Method.AIMNET2_WB97MD3]
NNP_METHODS = [Method.AIMNET2_WB97MD3]

XTBMethod = Literal[Method.GFN_FF, Method.GFN0_XTB, Method.GFN1_XTB, Method.GFN2_XTB]
XTB_METHODS = [Method.GFN_FF, Method.GFN0_XTB, Method.GFN1_XTB, Method.GFN2_XTB]

CompositeMethod = Literal[Method.HF3C, Method.B973C, Method.R2SCAN3C, Method.WB97X3C]
COMPOSITE_METHODS = [Method.HF3C, Method.B973C, Method.R2SCAN3C, Method.WB97X3C]

FFMethod = Literal[Method.OFF_SAGE_2_2_1]
FF_METHODS = [Method.OFF_SAGE_2_2_1]

PrepackagedMethod = XTBMethod | CompositeMethod | NNPMethod | FFMethod
PREPACKAGED_METHODS = [*XTB_METHODS, *COMPOSITE_METHODS, *NNP_METHODS, *FF_METHODS]

MethodWithCorrection = Literal[Method.WB97XD3, Method.WB97XV, Method.WB97MV, Method.WB97MD3BJ, Method.DSDBLYPD3BJ]
METHODS_WITH_CORRECTION = [Method.WB97XD3, Method.WB97XV, Method.WB97MV, Method.WB97MD3BJ, Method.DSDBLYPD3BJ, Method.B97D3BJ]
