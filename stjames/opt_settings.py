from typing import Sequence

from pydantic import PositiveFloat, PositiveInt

from .base import Base
from .constraint import Constraint


class OptimizationSettings(Base):
    """
    Geometry optimization settings.

    :param max_steps: maximum number of optimization steps
    :param transition_state: perform transition state optimization
    :param recalc_hess_every: recalculate Hessian every n steps (0 = never)
    :param max_gradient_threshold: convergence threshold for max gradient, in Hartree/Å
    :param rms_gradient_threshold: convergence threshold for RMS gradient, in Hartree/Å
    :param energy_threshold: convergence threshold for energy change, in Hartree
    :param optimize_cell: optimize unit cell (periodic systems only)
    :param constraints: geometric constraints to apply
    :param save_intermediate_steps: save intermediate geometries
    """

    max_steps: PositiveInt = 250
    transition_state: bool = False
    recalc_hess_every: int = 0

    max_gradient_threshold: PositiveFloat = 7e-4
    rms_gradient_threshold: PositiveFloat = 6e-4
    energy_threshold: PositiveFloat = 1e-6

    optimize_cell: bool = False

    constraints: Sequence[Constraint] = ()

    save_intermediate_steps: bool = True
