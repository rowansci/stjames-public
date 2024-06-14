from pydantic import PositiveFloat, PositiveInt

from .base import Base, LowercaseStrEnum


class ImplicitSolventModel(LowercaseStrEnum):
    CPCM = "cpcm"
    ALPB = "alpb"
    COSMO = "cosmo"
    NONE = "none"


class SolventSettings(Base):
    model: ImplicitSolventModel = ImplicitSolventModel.NONE
    epsilon: float = 78.36

    grid_points_per_atom: PositiveInt = 170
    vdw_scale: PositiveFloat = 1.2
    weight_cutoff: PositiveFloat = 1e-8
