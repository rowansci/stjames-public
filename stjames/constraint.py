from .base import Base, LowercaseStrEnum


class ConstraintType(LowercaseStrEnum):
    """Different sorts of constraints."""

    BOND = "bond"
    ANGLE = "angle"
    DIHEDRAL = "dihedral"


class Constraint(Base):
    """Represents a single constraint."""

    constraint_type: ConstraintType
    atoms: list[int]  # 1-indexed
