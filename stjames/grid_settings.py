import pydantic

from .base import Base, LowercaseStrEnum


class RadialGridType(LowercaseStrEnum):
    """What sort of radial grid (only one option for now)"""

    KRACK_KOSTER = "krack_koster"
    LMG = "lmg"


class GridSettings(Base):
    radial_grid_type: RadialGridType = RadialGridType.LMG
    angular_num_points: pydantic.PositiveInt = 434

    weight_cutoff: pydantic.PositiveFloat = 1e-11

    # for LMG
    radial_precision: pydantic.PositiveFloat = 1e-11

    # for other schemes, like KK
    radial_num_points: pydantic.PositiveInt = 75

    # pruning?
    prune: bool = True
    min_angular_points: pydantic.PositiveInt = 50
