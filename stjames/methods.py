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

    # hybrid GGA
    PBE0 = "pbe0"
    B3LYP = "b3lyp"
    B3PW91 = "b3pw91"
