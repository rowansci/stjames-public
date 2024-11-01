from typing import Optional, Self

from pydantic import PositiveFloat, PositiveInt, model_validator

from .base import Base, LowercaseStrEnum


class ConstraintType(LowercaseStrEnum):
    """Different sorts of constraints."""

    BOND = "bond"
    ANGLE = "angle"
    DIHEDRAL = "dihedral"
    FREEZE_ATOMS = "freeze_atoms"


class Constraint(Base):
    """
    Represents a single (absolute) constraint.

    :param constraint_type: which type
    :param atoms: the atoms in question. n.b. - these are 1-indexed!
    :param value: the value to constrain this to, leaving this blank sets the current value
    """

    constraint_type: ConstraintType
    atoms: list[PositiveInt]  # 1-indexed
    value: Optional[float] = None

    @model_validator(mode="after")
    def check_atom_list_length(self) -> Self:
        match self.constraint_type:
            case ConstraintType.BOND:
                if len(self.atoms) != 2:
                    raise ValueError("Bond constraint needs two atom indices!")
            case ConstraintType.ANGLE:
                if len(self.atoms) != 3:
                    raise ValueError("Angle constraint needs three atom indices!")
            case ConstraintType.DIHEDRAL:
                if len(self.atoms) != 4:
                    raise ValueError("Dihedral constraint needs four atom indices!")
            case ConstraintType.FREEZE_ATOMS:
                if len(self.atoms) == 0:
                    raise ValueError("Can't freeze atoms without any atoms to freeze!")
            case _:
                raise ValueError("Unknown constraint_type!")

        return self


class PairwiseHarmonicConstraint(Base):
    """
    Represents a harmonic constraint, with a characteristic spring constant.

    :param atoms: which atoms to apply to
    :param force_constant: the strength of the attraction, in kcal/mol/Å
    :param equilibrium: the distance at which force is zero
    """

    atoms: tuple[PositiveInt, PositiveInt]  # 1-indexed
    force_constant: PositiveFloat  # kcal/mol / Å**2
    equilibrium: PositiveFloat  # Å


class SphericalHarmonicConstraint(Base):
    """
    Represents a spherical harmonic constraint to keep a system near the origin.

    :param confining radius: the confining radius, in Å
    :param force_constant: the strength of the confinement, in kcal/mol/Å
    """

    confining_radius: PositiveFloat
    force_constant: PositiveFloat = 10  # kcal/mol / Å**2
