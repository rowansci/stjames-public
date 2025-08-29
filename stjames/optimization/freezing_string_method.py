"""Settings for the Freezing String Method (FSM)."""

from typing import Annotated, Self

from pydantic import AfterValidator, BaseModel, PositiveFloat, PositiveInt, model_validator

from ..base import LowercaseStrEnum, round_float


class FSMOptimizationCoordinates(LowercaseStrEnum):
    """Coordinate systems for FSM optimization step."""

    CARTESIAN = "cartesian"
    REDUNDANT_INTERNAL_COORDINATES = "redundant_internal_coordinates"


class FSMInterpolation(LowercaseStrEnum):
    """Methods for FSM interpolation step."""

    CARTESIAN = "cartesian"
    LINEAR_SYNCHRONOUS_TRANSIT = "linear_synchronous_transit"
    REDUNDANT_INTERNAL_COORDINATES = "redundant_internal_coordinates"


class FSMSettings(BaseModel):
    """
    Settings for the Freezing String Method (FSM) TS search.

    :param opt_coords: coordinate system to use for optimization
    :param interpolation_method: method to use for interpolation between nodes
    :param min_num_nodes: minimum number of nodes to use in the string
    :param num_interpolation_points: number of interpolation points to use between end nodes
    :param max_optimizer_iterations: maximum number of optimizer iterations to perform  (scipy.minimize maxiter)
    :param max_line_search_steps: maximum number of line search steps to perform (scipy.minimize maxls)
    :param max_displacement: maximum displacement for a single coordinate (Å for distances, internally converted for angles)
    """

    optimization_coordinates: FSMOptimizationCoordinates = FSMOptimizationCoordinates.CARTESIAN
    interpolation_method: FSMInterpolation = FSMInterpolation.REDUNDANT_INTERNAL_COORDINATES

    min_num_nodes: PositiveInt = 18
    num_interpolation_points: PositiveInt = 100

    max_optimizer_iterations: PositiveInt = 3
    max_line_search_steps: PositiveInt = 2
    max_displacement: Annotated[PositiveFloat, AfterValidator(round_float(3))] = 0.3

    @model_validator(mode="after")
    def sanity_check(self) -> Self:
        """Sanity check for settings."""
        if self.min_num_nodes < 5:
            raise ValueError("Must have at least 5 nodes for the Freezing String Method.")
        if self.max_displacement < 0.001:
            raise ValueError("Maximum displacement must be at least 0.001 Å.")

        return self
