from typing import Optional
import pydantic

from .base import Base, LowercaseStrEnum
from .int_settings import IntSettings
from .grid_settings import GridSettings
from .diis_settings import DIISSettings


class SCFInitMethod(LowercaseStrEnum):
    # (See https://manual.q-chem.com/5.2/Ch4.S4.SS2.html for a nice overview here.)
    SAD = "sad"
    CORE = "core"
    #    GWH = "gwh"
    READ = "read"


class OrthonormalizationMethod(LowercaseStrEnum):
    SYMMETRIC = "symmetric"
    CANONICAL = "canonical"
    # cholesky, in future?


class SCFSettings(Base):
    max_iters: int = 100
    init_method: SCFInitMethod = SCFInitMethod.SAD

    int_settings: IntSettings = IntSettings()
    grid_settings: GridSettings = GridSettings()
    diis_settings: DIISSettings = DIISSettings()

    #### orthonormalization
    orthonormalization: OrthonormalizationMethod = OrthonormalizationMethod.CANONICAL

    #### damping
    do_damping: bool = True
    # when should we stop damping?
    end_damping_error: pydantic.PositiveFloat = 0.5
    # what damping factor should we use?
    damping_factor: float = pydantic.Field(ge=0, le=1, default=0.5)

    #### level shifting
    do_level_shift: bool = True
    # how much? (Eh)
    level_shift_magnitude: pydantic.PositiveFloat = 0.1
    # when should we stop?
    end_level_shift_error: pydantic.PositiveFloat = 0.1

    #### reset incremental fock build
    rebuild_frequency: pydantic.PositiveInt = 20

    #### when are we converged?
    energy_threshold: pydantic.PositiveFloat = 1e-6
    rms_error_threshold: pydantic.PositiveFloat = 1e-8
    max_error_threshold: pydantic.PositiveFloat = 1e-5

    #### DIIS
    do_diis: bool = True
    # error below which we'll start DIIS
    start_diis_max_error: pydantic.PositiveFloat = 0.5
    # first iteration we'll consider starting DIIS
    start_diis_iter: pydantic.PositiveInt = 3
    # iteration past which we'll start DIIS even if error is high
    start_diis_anyway: pydantic.PositiveInt = 7

    # if ``read`` initialization is selected
    initial_density_matrix_guess: Optional[list[list[float]]] = None
