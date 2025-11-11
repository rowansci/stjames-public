import re
from typing import Annotated, Self, TypeAlias

import numpy as np
import pydantic

from .base import Base
from .types import Matrix3x3, round_matrix3x3

Bool3: TypeAlias = tuple[bool, bool, bool]


class PeriodicCell(Base):
    lattice_vectors: Annotated[Matrix3x3, pydantic.AfterValidator(round_matrix3x3(6))]
    is_periodic: Bool3 = (True, True, True)

    def __repr__(self) -> str:
        """
        Return a string representation of the PeriodicCell.

        >>> PeriodicCell(lattice_vectors=[[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        PeriodicCell(lattice_vectors=((1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)), is_periodic=(True, True, True))
        """
        return f"PeriodicCell(lattice_vectors={self.lattice_vectors}, is_periodic={self.is_periodic})"

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

    @classmethod
    def from_string(cls, string: str) -> Self:
        """
        Create a PeriodicCell from a string representation.

        >>> PeriodicCell.from_string("[[1, 2, 3], [4, 5, 6], [7, 8, 9]]")
        PeriodicCell(lattice_vectors=((1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)), is_periodic=(True, True, True))
        >>> PeriodicCell.from_string("[(1, -2, 3.0), [4e1, +5, 6], (7.1, 8, 9)]")
        PeriodicCell(lattice_vectors=((1.0, -2.0, 3.0), (40.0, 5.0, 6.0), (7.1, 8.0, 9.0)), is_periodic=(True, True, True))
        """
        if match := re.match(_cell_pattern, string.strip()):
            cell = [[float(num) for num in line.split(",")] for line in match.groups()]
        else:
            raise ValueError(f"String '{string}' is not a valid representation of a PeriodicCell.")

        return cls(lattice_vectors=cell)


_ob = r"\s*[\[\(]\s*"
_cb = r"\s*[\]\)]\s*"
_line = rf"\s*{_ob}(.+?){_cb}\s*"
_cell_pattern = f"{_ob}{_line},{_line},{_line}{_cb}"
