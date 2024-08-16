from typing import Optional, Self

from pydantic import PositiveFloat, PositiveInt, model_validator

from .base import Base


class BasisSetOverride(Base):
    name: str
    atomic_numbers: Optional[list[PositiveInt]] = None
    atoms: Optional[list[PositiveInt]] = None  # 1-indexed

    @model_validator(mode="after")
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
    cutoff_threshold: PositiveFloat = 1e-10
