from pydantic import PositiveFloat, PositiveInt

from .base import Base, LowercaseStrEnum


class ConstraintType(LowercaseStrEnum):
    """Different sorts of constraints."""

    BOND = "bond"
    ANGLE = "angle"
    DIHEDRAL = "dihedral"


class Constraint(Base):
    """Represents a single (absolute) constraint."""

    constraint_type: ConstraintType
    atoms: list[PositiveInt]  # 1-indexed


class PairwiseHarmonicConstraint(Base):
    """
    Represents a harmonic constraint, with a characteristic spring constant.

    :param atoms: whch atoms to apply to
    :force_constant: the strength of the attraction
    :equilibrium: the distance at which force is zero
    """

    atoms: tuple[PositiveInt, PositiveInt]  # 1-indexed
    force_constant: PositiveFloat  # kcal/mol / Å**2
    equilibrium: PositiveFloat  # Å


class SphericalHarmonicConstraint(Base):
    """
    Represents a spherical harmonic constraint to keep a system near the origin.
    """

    confining_radius: PositiveFloat
    force_constant: PositiveFloat = 10  # kcal/mol / Å**2
