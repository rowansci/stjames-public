from typing import Self

from pydantic import PositiveFloat, PositiveInt, model_validator

from .base import Base


class BasisSetOverride(Base):
    """
    Basis set override for specific atoms or elements.

    :param name: basis set name for override
    :param atomic_numbers: atomic numbers to override (mutually exclusive with atoms)
    :param atoms: 1-indexed atom indices to override (mutually exclusive with atomic_numbers)
    """

    name: str
    atomic_numbers: list[PositiveInt] | None = None
    atoms: list[PositiveInt] | None = None

    @model_validator(mode="after")
    def check_override(self) -> Self:
        # ^ is xor
        assert (self.atomic_numbers is not None) ^ (self.atoms is not None), "Exactly one of atomic_numbers or atoms must be specified!"
        return self


class BasisSet(Base):
    """
    Atomic orbital basis set.

    :param name: basis set name (e.g., def2-SVP, cc-pVDZ)
    :param overrides: element or atom-specific basis set overrides
    :param cutoff_threshold: basis function screening threshold
    """

    name: str

    # do we want to override the default basis set for specific atoms or elements?
    overrides: list[BasisSetOverride] | None = []

    # value below which a basis function can be ignored
    # (for improving DFT grid calcs, as per Stratmann/Scuseria/Frisch CPL 1996)
    # this shouldn't really need to be modified...
    cutoff_threshold: PositiveFloat = 1e-10
