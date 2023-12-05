import pydantic

from .base import Base, LowercaseStrEnum


class RadialGridType(LowercaseStrEnum):
    """What sort of radial grid (only one option for now)"""

    KRACK_KOSTER = "krack_koster"
    LMG = "lmg"


class GridSettings(Base):
    radial_grid_type: RadialGridType = RadialGridType.LMG
    angular_num_points: pydantic.PositiveInt = 590

    weight_cutoff: pydantic.PositiveFloat = 1e-14

    # for LMG
    radial_precision: pydantic.PositiveFloat = 1e-12

    # for other schemes, like KK
    radial_num_points: pydantic.PositiveInt = 99

    # pruning?
    prune: bool = True
    min_angular_points: pydantic.PositiveInt = 50
