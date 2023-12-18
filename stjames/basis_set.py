import pydantic
from typing import Optional, Self

from .base import Base


class BasisSetOverride(Base):
    name: str
    atomic_numbers: Optional[list[pydantic.PositiveInt]] = None
    atoms: Optional[list[pydantic.PositiveInt]] = None  # 1-indexed

    @pydantic.model_validator(mode="after")
    def check_override(self) -> Self:
        # ^ is xor
        assert (self.atomic_numbers is not None) ^ (self.atoms is not None), "Exactly one of ``atomic_numbers`` or ``atoms`` must be specified!"
        return self


class BasisSet(Base):
    name: str

    # do we want to override the default basis set for specific atoms or elements?
    overrides: Optional[list[BasisSetOverride]] = []

    # value below which a basis function can be ignored
    # (for improving DFT grid calcs, as per Stratmann/Scuseria/Frisch CPL 1996)
    # this shouldn't really need to be modified...
    cutoff_threshold: pydantic.PositiveFloat = 1e-10
