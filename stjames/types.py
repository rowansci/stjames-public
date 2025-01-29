from typing import Callable, TypeAlias

UUID: TypeAlias = str

Vector3D: TypeAlias = tuple[float, float, float]
Vector3DPerAtom: TypeAlias = list[Vector3D]

FloatPerAtom: TypeAlias = list[float]

Matrix3x3: TypeAlias = tuple[Vector3D, Vector3D, Vector3D]


def round_vector3d(precision: int = 3) -> Callable[[Vector3D], Vector3D]:
    """Create a validator that rounds each component of a Vector3D to the specified precision."""

    def rounder(vector: Vector3D) -> Vector3D:
        return (round(vector[0], precision), round(vector[1], precision), round(vector[2], precision))

    return rounder


def round_vector3d_per_atom(precision: int = 3) -> Callable[[Vector3DPerAtom], Vector3DPerAtom]:
    """Create a validator that rounds each vector in Vector3DPerAtom to the specified precision."""
    vector_rounder = round_vector3d(precision)

    def rounder(vectors: Vector3DPerAtom) -> Vector3DPerAtom:
        return [vector_rounder(v) for v in vectors]

    return rounder
