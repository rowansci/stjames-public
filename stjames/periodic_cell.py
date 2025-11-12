import ast
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

        >>> PeriodicCell.from_string("[(1, -2, 3.0), [4e1, +5, 6], (7.1, 8, 9)]")
        PeriodicCell(lattice_vectors=((1.0, -2.0, 3.0), (40.0, 5.0, 6.0), (7.1, 8.0, 9.0)), is_periodic=(True, True, True))
        >>> PeriodicCell.from_string("[(a, -2, 3.0), [4e1, +5, 6], (7.1, 8, 9)]")
        Traceback (most recent call last):
        ...
        ValueError: Could not parse cell from string: [(a, -2, 3.0), [4e1, +5, 6], (7.1, 8, 9)]
        >>> PeriodicCell.from_string("[(1, -2, 3.0), [4e1, +5, 6], (7.1, 8, 9, 10)]")
        Traceback (most recent call last):
        ...
        ValueError: Cell must be a 3x3 matrix: got [(1, -2, 3.0), [4e1, +5, 6], (7.1, 8, 9, 10)]
        >>> PeriodicCell.from_string("[('a', -2, 3.0), [4e1, +5, 6], (7.1, 8, 9)]")
        Traceback (most recent call last):
        ValueError: Cell must be a 3x3 matrix of numbers: got [('a', -2, 3.0), [4e1, +5, 6], (7.1, 8, 9)]

        """
        try:
            cell = ast.literal_eval(string)
        except (ValueError, SyntaxError) as e:
            raise ValueError(f"Could not parse cell from string: {string}") from e

        if not isinstance(cell, (list, tuple)):
            raise ValueError(f"Cell must be a list or tuple: got {string}")
        if not len(cell) == 3:
            raise ValueError(f"Cell must be a 3x3 matrix: got {string}")
        if not all(isinstance(row, (list, tuple)) for row in cell):
            raise ValueError(f"Cell must be a 3x3 matrix: got {string}")
        if not all(len(row) == 3 for row in cell):
            raise ValueError(f"Cell must be a 3x3 matrix: got {string}")
        if not all(all(isinstance(x, (int, float)) for x in row) for row in cell):
            raise ValueError(f"Cell must be a 3x3 matrix of numbers: got {string}")

        return cls(lattice_vectors=cell)
