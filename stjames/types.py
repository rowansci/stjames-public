from typing import TypeAlias

UUID: TypeAlias = str

Vector3D: TypeAlias = tuple[float, float, float]
Vector3DPerAtom: TypeAlias = list[Vector3D]

Matrix3x3: TypeAlias = tuple[Vector3D, Vector3D, Vector3D]
