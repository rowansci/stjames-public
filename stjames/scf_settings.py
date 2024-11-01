from typing import Any, Optional

import pydantic
from pydantic import PositiveFloat, PositiveInt

from .base import Base, LowercaseStrEnum
from .diis_settings import DIISSettings
from .grid_settings import GridSettings
from .int_settings import IntSettings


class SCFInitMethod(LowercaseStrEnum):
    # (See https://manual.q-chem.com/5.2/Ch4.S4.SS2.html for a nice overview here.)
    SAD = "sad"
    CORE = "core"
    #    GWH = "gwh"
    READ = "read"


class OrthonormalizationMethod(LowercaseStrEnum):
    SYMMETRIC = "symmetric"
    CANONICAL = "canonical"
    # Cholesky, in future?


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
    end_damping_error: PositiveFloat = 0.1
    # what damping factor should we use?
    damping_factor: float = pydantic.Field(ge=0, le=1, default=0.7)

    #### level shifting
    do_level_shift: bool = True
    # how much? (Eh)
    level_shift_magnitude: PositiveFloat = 0.25
    # when should we stop?
    end_level_shift_error: PositiveFloat = 0.1

    #### incremental
    # do incremental fock build?
    do_incremental: bool = True
    # reset incremental fock build
    rebuild_frequency: PositiveInt = 20

    #### when are we converged?
    energy_threshold: PositiveFloat = 1e-6
    rms_error_threshold: PositiveFloat = 1e-8
    max_error_threshold: PositiveFloat = 1e-5

    #### DIIS
    do_diis: bool = True
    # error below which we'll start DIIS
    start_diis_max_error: pydantic.PositiveFloat = 0.2
    # first iteration we'll consider starting DIIS
    start_diis_iter: PositiveInt = 3
    # iteration past which we'll start DIIS even if error is high
    start_diis_anyway: PositiveInt = 7

    # if ``read`` initialization is selected
    initial_density_matrix_guess: Optional[list[list[float]]] = None

    def model_post_init(self, __context: Any) -> None:
        # disable incremental Fock for RI
        if self.int_settings.resolution_of_the_identity:
            self.do_incremental = False
