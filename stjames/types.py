from typing import Callable, Iterable, Protocol, TypeAlias, overload

UUID: TypeAlias = str

Vector3D: TypeAlias = tuple[float, float, float]
Vector3DPerAtom: TypeAlias = list[Vector3D]

FloatPerAtom: TypeAlias = list[float]

Matrix3x3: TypeAlias = tuple[Vector3D, Vector3D, Vector3D]


class _Iterable_Rounder(Protocol):
    """Callable with return types conditioned on the input type(s)."""

    @overload
    def __call__(self, values: None, /) -> None: ...
    @overload
    def __call__(self, values: Iterable[float], /) -> list[float]: ...
    @overload
    def __call__(self, values: Iterable[float | None], /) -> list[float | None]: ...


def round_list(round_to: int = 6) -> _Iterable_Rounder:
    """Round values (even if Nones), preserving types."""

    @overload
    def rounder(values: None, /) -> None: ...
    @overload
    def rounder(values: Iterable[float], /) -> list[float]: ...
    @overload
    def rounder(values: Iterable[float | None], /) -> list[float | None]: ...

    def rounder(values: Iterable[float] | Iterable[float | None] | None, /) -> list[float] | list[float | None] | None:
        if values is None:
            return None

        return [round(v, round_to) if v is not None else None for v in values]

    return rounder


def round_vector3d(round_to: int = 6) -> Callable[[Vector3D], Vector3D]:
    """Create a validator that rounds each component of a Vector3D to a given number of decimal places."""

    def rounder(vector: Vector3D) -> Vector3D:
        return (round(vector[0], round_to), round(vector[1], round_to), round(vector[2], round_to))

    return rounder


def round_optional_vector3d(round_to: int = 6) -> Callable[[Vector3D | None], Vector3D | None]:
    """Create a validator that rounds each component of a Vector3D to a given number of decimal places, handling None."""

    def rounder(vector: Vector3D | None) -> Vector3D | None:
        if vector is None:
            return None
        return (round(vector[0], round_to), round(vector[1], round_to), round(vector[2], round_to))

    return rounder


def round_vector3d_per_atom(round_to: int = 6) -> Callable[[Vector3DPerAtom], Vector3DPerAtom]:
    """Create a validator that rounds each vector in Vector3DPerAtom to a given number of decimal places."""
    vector_rounder = round_vector3d(round_to)

    def rounder(vectors: Vector3DPerAtom) -> Vector3DPerAtom:
        return [vector_rounder(v) for v in vectors]

    return rounder


def round_optional_vector3d_per_atom(round_to: int = 6) -> Callable[[Vector3DPerAtom | None], Vector3DPerAtom | None]:
    """Create a validator that rounds each vector in Vector3DPerAtom to a given number of decimal places, handling None."""
    vector_rounder = round_vector3d(round_to)

    def rounder(vectors: Vector3DPerAtom | None) -> Vector3DPerAtom | None:
        if vectors is None:
            return None
        return [vector_rounder(v) for v in vectors]

    return rounder


def round_matrix3x3(round_to: int = 6) -> Callable[[Matrix3x3], Matrix3x3]:
    """Create a validator that rounds each vector in a Matrix3x3 to a given number of decimal places."""

    # Use the round_vector3d function to round each Vector3D in the Matrix3x3
    vector_rounder = round_vector3d(round_to)

    def rounder(matrix: Matrix3x3) -> Matrix3x3:
        return (vector_rounder(matrix[0]), vector_rounder(matrix[1]), vector_rounder(matrix[2]))

    return rounder


def round_optional_matrix3x3(round_to: int = 3) -> Callable[[Matrix3x3 | None], Matrix3x3 | None]:
    """Create a validator that rounds each vector in an Optional Matrix3x3 to a given number of decimal places."""

    # Use the round_vector3d function to round each Vector3D in the Matrix3x3
    vector_rounder = round_vector3d(round_to)

    def rounder(matrix: Matrix3x3 | None) -> Matrix3x3 | None:
        if matrix is None:
            return None
        return (vector_rounder(matrix[0]), vector_rounder(matrix[1]), vector_rounder(matrix[2]))

    return rounder


def round_optional_float_per_atom(round_to: int = 6) -> Callable[[FloatPerAtom | None], FloatPerAtom | None]:
    """Create a validator that rounds each float in FloatPerAtom to a given number of decimal places, handling None."""

    def rounder(values: FloatPerAtom | None) -> FloatPerAtom | None:
        if values is None:
            return None
        return [round(value, round_to) for value in values]

    return rounder


def round_float_per_atom(round_to: int = 6) -> Callable[[FloatPerAtom], FloatPerAtom]:
    """Create a validator that rounds each float in FloatPerAtom to a given number of decimal places, handling None."""

    def rounder(values: FloatPerAtom) -> FloatPerAtom:
        return [round(value, round_to) for value in values]

    return rounder
