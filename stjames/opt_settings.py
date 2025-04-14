from typing import Sequence

from pydantic import PositiveFloat, PositiveInt

from .base import Base
from .constraint import Constraint


class OptimizationSettings(Base):
    max_steps: PositiveInt = 250
    transition_state: bool = False
    recalc_hess_every: int = 0  # Recalculate the hessian every `n` steps (or never with 0)

    # when are we converged? (Hartree and Hartree/Å)
    max_gradient_threshold: PositiveFloat = 7e-4
    rms_gradient_threshold: PositiveFloat = 6e-4
    energy_threshold: PositiveFloat = 1e-6

    # for periodic systems only
    optimize_cell: bool = False

    constraints: Sequence[Constraint] = tuple()
