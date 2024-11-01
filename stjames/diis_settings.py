import pydantic

from .base import Base, LowercaseStrEnum


class DIISStrategy(LowercaseStrEnum):
    # regular Pulay DIIS
    DIIS = "diis"

    # Hu/Yang JCP 2010 - ADIIS
    ADIIS = "adiis"

    # first ADIIS, then DIIS
    ADIIS_DIIS = "adiis_diis"


class DIISSettings(Base):
    strategy: DIISStrategy = DIISStrategy.ADIIS_DIIS
    subspace_size: pydantic.PositiveInt = 12

    # if it's a hybrid strategy, where do we transition?
    adiis_diis_blend_start: pydantic.PositiveFloat = 1e-1
    adiis_diis_blend_stop: pydantic.PositiveFloat = 1e-4
