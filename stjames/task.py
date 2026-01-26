from .base import LowercaseStrEnum


class Task(LowercaseStrEnum):
    """Calculation task type."""

    ENERGY = "energy"
    GRADIENT = "gradient"
    OPTIMIZE = "optimize"
    OPTIMIZE_TS = "optimize_ts"
    CHARGE = "charge"
    SPIN_DENSITY = "spin_density"
    DIPOLE = "dipole"
    HESSIAN = "hessian"
    FREQUENCIES = "frequencies"
    STRESS = "stress"
