from pydantic import Field, PositiveFloat, PositiveInt

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


class OptimizationSettings(Base):
    max_steps: PositiveInt = 250
    transition_state: bool = False

    # when are we converged?
    max_gradient_threshold: PositiveFloat = 4.5e-4
    rms_gradient_threshold: PositiveFloat = 3.0e-4
    energy_threshold: PositiveFloat = 1e-6

    constraints: list[Constraint] = Field(default_factory=list)
