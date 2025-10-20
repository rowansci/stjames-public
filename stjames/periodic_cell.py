from typing import Annotated, TypeAlias

import numpy as np
import pydantic

from .base import Base
from .types import Matrix3x3, round_matrix3x3

Bool3: TypeAlias = tuple[bool, bool, bool]


class PeriodicCell(Base):
    lattice_vectors: Annotated[Matrix3x3, pydantic.AfterValidator(round_matrix3x3(6))]
    is_periodic: Bool3 = (True, True, True)

    @pydantic.field_validator("lattice_vectors")
    @classmethod
    def check_tensor_3D(cls, v: Matrix3x3) -> Matrix3x3:
        if len(v) != 3 or any(len(row) != 3 for row in v):
            raise ValueError("Cell tensor must be a 3x3 list of floats")

        return v

    @pydantic.field_validator("is_periodic")
    @classmethod
    def check_pbc(cls, v: Bool3) -> Bool3:
        if not any(v):
            raise ValueError("For periodic boundary conditions, at least one dimension must be periodic!")
        return v

    @pydantic.computed_field  # type: ignore[misc, prop-decorator, unused-ignore]
    @property
    def volume(self) -> float:
        return float(np.abs(np.linalg.det(np.array(self.lattice_vectors))))
