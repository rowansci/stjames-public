from pytest import raises

from stjames import Constraint


def test_invalid_constraint_settings() -> None:
    Constraint(constraint_type="bond", atoms=(1, 2))

    with raises(ValueError):
        Constraint(constraint_type="bond", atoms=(1, 2, 3))

    with raises(ValueError):
        Constraint(constraint_type="bond", atoms=(1, 2, 3, 4))

    Constraint(constraint_type="angle", atoms=(1, 2, 3))

    with raises(ValueError):
        Constraint(constraint_type="angle", atoms=(1, 2))

    with raises(ValueError):
        Constraint(constraint_type="angle", atoms=(1, 2, 3, 4))

    Constraint(constraint_type="dihedral", atoms=(1, 2, 3, 4))

    with raises(ValueError):
        Constraint(constraint_type="dihedral", atoms=(1, 2))

    with raises(ValueError):
        Constraint(constraint_type="dihedral", atoms=(1, 2, 3))
