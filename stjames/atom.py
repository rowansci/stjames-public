from typing import Annotated, Self, Sequence

from pydantic import AfterValidator, NonNegativeInt

from .base import Base
from .data import ELEMENT_SYMBOL, SYMBOL_ELEMENT
from .types import Vector3D, round_vector3d


class Atom(Base):
    atomic_number: NonNegativeInt
    position: Annotated[Vector3D, AfterValidator(round_vector3d(8))]  # in Å

    def __repr__(self) -> str:
        """
        >>> Atom(atomic_number=2, position=[0, 1, 2])
        Atom(2, [0.00000, 1.00000, 2.00000])
        """
        x, y, z = self.position
        return f"Atom({self.atomic_number}, [{x:.5f}, {y:.5f}, {z:.5f}])"

    def __str__(self) -> str:
        """
        >>> str(Atom(atomic_number=2, position=[0, 1, 2]))
        'He    0.0000000000    1.0000000000    2.0000000000'
        """
        x, y, z = self.position
        return f"{self.atomic_symbol:2} {x:15.10f} {y:15.10f} {z:15.10f}"

    @property
    def atomic_symbol(self) -> str:
        """
        >>> Atom(atomic_number=2, position=[0, 1, 2]).atomic_symbol
        'He'
        """
        return ELEMENT_SYMBOL[self.atomic_number]

    def edited(self, atomic_number: int | None = None, position: Sequence[float] | None = None) -> Self:
        """
        Create a new Atom with the specified changes.

        >>> a = Atom(atomic_number=2, position=[0, 1, 2])
        >>> a2 = a.edited(3)
        >>> a is a2
        False
        >>> a2
        Atom(3, [0.00000, 1.00000, 2.00000])
        """
        if atomic_number is None:
            atomic_number = self.atomic_number
        if position is None:
            position = list(self.position)

        return self.__class__(atomic_number=atomic_number, position=position)

    @classmethod
    def from_xyz(cls: type[Self], xyz_line: str) -> Self:
        """
        >>> Atom.from_xyz("H 0 0 0")
        Atom(1, [0.00000, 0.00000, 0.00000])
        """
        name, *xyz = xyz_line.split()
        symbol = int(name) if name.isdigit() else SYMBOL_ELEMENT[name.title()]
        if not len(xyz) == 3:
            raise ValueError("XYZ file should have 3 coordinates per atom")
        return cls(atomic_number=symbol, position=xyz)
