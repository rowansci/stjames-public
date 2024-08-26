from pydantic import PositiveFloat

from .base import Base, LowercaseStrEnum


class ConstraintType(LowercaseStrEnum):
    """Different sorts of constraints."""

    BOND = "bond"
    ANGLE = "angle"
    DIHEDRAL = "dihedral"


class Constraint(Base):
    """Represents a single (absolute) constraint."""

    constraint_type: ConstraintType
    atoms: list[int]  # 1-indexed


class PairwiseHarmonicConstraint(Base):
    """
    Represents a harmonic constraint, with a characteristic spring constant.
    """

    atoms: tuple[int, int]  # 1-indexed
    spring_constant: PositiveFloat  # kcal/mol / Å**2
