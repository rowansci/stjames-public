import pydantic
from typing import Self

from .base import Base, LowercaseStrEnum


class ConstraintType(LowercaseStrEnum):
    """Different sorts of constraints."""

    BOND = "bond"
    ANGLE = "angle"
    DIHEDRAL = "dihedral"


class Constraint(Base):
    """Represents a single constraint."""

    constraint_type: ConstraintType
    atoms: list[int]


class OptimizationAlgorithm(LowercaseStrEnum):
    """Which algorithm to use."""

    # default, L-BGFS
    L_BGFS = "l_bgfs"

    # P-RFO, for TS
    P_RFO = "p_rfo"

    # Polak-Ribiere
    CONJUGATE_GRADIENT = "conjugate_gradient"

    # steepest descent
    STEEPEST_DESCENT = "steepest_descent"


class OptimizationSettings(Base):
    max_steps: pydantic.PositiveInt = 100
    algorithm: OptimizationAlgorithm = OptimizationAlgorithm.L_BGFS
    max_stepsize: pydantic.PositiveFloat = 0.05

    # when are we converged?
    max_gradient_threshold: pydantic.PositiveFloat = 4.5e-4
    rms_gradient_threshold: pydantic.PositiveFloat = 3.0e-4
    energy_threshold: pydantic.PositiveFloat = 1e-6

    constraints: list[Constraint] = pydantic.Field(default_factory=list)

    @pydantic.model_validator(mode="after")
    def check_rms_max_grad(self) -> Self:
        rms_grad = self.rms_gradient_threshold
        max_grad = self.max_gradient_threshold
        assert abs(max_grad / rms_grad - 1.50) < 1e-3, "DLFIND hard-codes this relationship between RMS and max gradient!"
        return self
