import pydantic

from .base import Base, LowercaseStrEnum


class ImplicitSolventModel(LowercaseStrEnum):
    CPCM = "cpcm"
    ALPB = "alpb"
    COSMO = "cosmo"
    NONE = "none"


class SolventSettings(Base):
    model: ImplicitSolventModel = ImplicitSolventModel.NONE
    epsilon: float = 78.36

    grid_points_per_atom: pydantic.PositiveInt = 170
    vdw_scale: pydantic.PositiveFloat = 1.2
    weight_cutoff: pydantic.PositiveFloat = 1e-8
