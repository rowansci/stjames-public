from typing import Literal

from .base import LowercaseStrEnum
from .engine import Engine


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
    MACE_MP_0B2_L = "mace_mp_0b2_l"
    OMOL25_CONSERVING_S = "omol25_conserving_s"
    UMA_M_OMOL = "uma_m_omol"
    UMA_S_OMOL = "uma_s_omol"
    ORB_V3_CONSERVATIVE_INF_OMAT = "orb_v3_conservative_inf_omat"
    ORB_V3_CONSERVATIVE_OMOL = "orb_v3_conservative_omol"

    GFN_FF = "gfn_ff"
    GFN0_XTB = "gfn0_xtb"
    GFN1_XTB = "gfn1_xtb"
    GFN2_XTB = "gfn2_xtb"
    G_XTB = "g_xtb"

    # this was going to be removed, but Jonathon wrote such a nice basis set test... it's off the front end.
    BP86 = "bp86"

    OFF_SAGE_2_2_1 = "off_sage_2_2_1"

    EGRET_1 = "egret_1"
    EGRET_1E = "egret_1e"
    EGRET_1T = "egret_1t"

    def default_engine(self, *, is_periodic: bool = False) -> Engine:
        """
        Return the canonical Engine for this quantum-chemistry method.

        :param bool is_periodic:
            If ``True`` **and** the method is in the XTB family, return
            ``"tblite"`` (periodic-capable backend) instead of ``"xtb"``.
        :returns: lower-case engine identifier (e.g. ``"psi4"``, ``"mace"``).

        >>> Method.MACE_MP_0B2_L.default_engine().value
        'mace'
        >>> Method.GFN2_XTB.default_engine().value
        'xtb'
        >>> Method.GFN2_XTB.default_engine(is_periodic=True).value
        'tblite'
        """
        match self:
            case Method.AIMNET2_WB97MD3:
                return Engine.AIMNET2
            case Method.MACE_MP_0B2_L:
                return Engine.MACE
            case Method.OMOL25_CONSERVING_S | Method.UMA_M_OMOL | Method.UMA_S_OMOL:
                return Engine.OMOL25
            case Method.ORB_V3_CONSERVATIVE_INF_OMAT | Method.ORB_V3_CONSERVATIVE_OMOL:
                return Engine.ORB
            case method if method in XTB_METHODS:
                return Engine.TBLITE if is_periodic else Engine.XTB
            case Method.OFF_SAGE_2_2_1:
                return Engine.OPENFF
            case Method.EGRET_1 | Method.EGRET_1E | Method.EGRET_1T:
                return Engine.EGRET
            case _:
                # All remaining methods (HF, DFT, composite, etc.) fall back to Psi4
                return Engine.PSI4


PrepackagedNNPMethod = Literal[
    Method.AIMNET2_WB97MD3,
    Method.OMOL25_CONSERVING_S,
    Method.UMA_M_OMOL,
    Method.UMA_S_OMOL,
    Method.ORB_V3_CONSERVATIVE_INF_OMAT,
    Method.ORB_V3_CONSERVATIVE_OMOL,
    Method.EGRET_1,
    Method.EGRET_1E,
    Method.EGRET_1T,
]

PREPACKAGED_NNP_METHODS = [
    Method.AIMNET2_WB97MD3,
    Method.OMOL25_CONSERVING_S,
    Method.UMA_M_OMOL,
    Method.UMA_S_OMOL,
    Method.ORB_V3_CONSERVATIVE_INF_OMAT,
    Method.ORB_V3_CONSERVATIVE_OMOL,
    Method.EGRET_1,
    Method.EGRET_1E,
    Method.EGRET_1T,
]

CorrectableNNPMethod = Literal[Method.MACE_MP_0B2_L]
CORRECTABLE_NNP_METHODS = [Method.MACE_MP_0B2_L]

NNPMethod = PrepackagedNNPMethod | CorrectableNNPMethod
NNP_METHODS = [*PREPACKAGED_NNP_METHODS, *CORRECTABLE_NNP_METHODS]

XTBMethod = Literal[Method.GFN_FF, Method.GFN0_XTB, Method.GFN1_XTB, Method.GFN2_XTB, Method.G_XTB]
XTB_METHODS = [Method.GFN_FF, Method.GFN0_XTB, Method.GFN1_XTB, Method.GFN2_XTB, Method.G_XTB]

CompositeMethod = Literal[Method.HF3C, Method.B973C, Method.R2SCAN3C, Method.WB97X3C]
COMPOSITE_METHODS = [Method.HF3C, Method.B973C, Method.R2SCAN3C, Method.WB97X3C]

FFMethod = Literal[Method.OFF_SAGE_2_2_1]
FF_METHODS = [Method.OFF_SAGE_2_2_1]

PrepackagedMethod = XTBMethod | CompositeMethod | PrepackagedNNPMethod | FFMethod
PREPACKAGED_METHODS = [*XTB_METHODS, *COMPOSITE_METHODS, *PREPACKAGED_NNP_METHODS, *FF_METHODS]

MethodWithCorrection = Literal[Method.WB97XD3, Method.WB97XV, Method.WB97MV, Method.WB97MD3BJ, Method.DSDBLYPD3BJ]
METHODS_WITH_CORRECTION = [Method.WB97XD3, Method.WB97XV, Method.WB97MV, Method.WB97MD3BJ, Method.DSDBLYPD3BJ, Method.B97D3BJ]

MGGAFunctionals = Literal[Method.R2SCAN, Method.R2SCAN3C, Method.TPSS, Method.TPSSH, Method.M06L, Method.M06, Method.M062X, Method.WB97MD3BJ, Method.WB97MV]
MGGA_FUNCTIONALS = [Method.R2SCAN, Method.R2SCAN3C, Method.TPSS, Method.TPSSH, Method.M06L, Method.M06, Method.M062X, Method.WB97MD3BJ, Method.WB97MV]
