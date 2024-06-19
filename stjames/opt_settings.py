from pydantic import Field, PositiveFloat, PositiveInt

from .base import Base
from .constraint import Constraint


class OptimizationSettings(Base):
    max_steps: PositiveInt = 250
    transition_state: bool = False

    # when are we converged?
    max_gradient_threshold: PositiveFloat = 4.5e-4
    rms_gradient_threshold: PositiveFloat = 3.0e-4
    energy_threshold: PositiveFloat = 1e-6

    constraints: list[Constraint] = Field(default_factory=list)
