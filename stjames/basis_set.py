import pydantic

from .base import Base


class BasisSet(Base):
    name: str

    # value below which a basis function can be ignored
    # (for improving DFT grid calcs, as per Stratmann/Scuseria/Frisch CPL 1996)
    # this shouldn't really need to be modified...
    cutoff_threshold: pydantic.PositiveFloat = 1e-10
