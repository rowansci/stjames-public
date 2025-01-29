from stjames import VibrationalMode


def test_vector_rounding() -> None:
    test_vectors = [(1.23456789345346, 2.34567890346346, 3.45678901346346436), (4.567890126343636, 5.678901236346436, 6.7890123436436436)]
    mode = VibrationalMode(frequency=1.0, reduced_mass=1.0, force_constant=1.0, displacements=test_vectors)

    expected_vectors = [(1.234568, 2.345679, 3.456789), (4.567890, 5.678901, 6.789012)]
    assert mode.displacements == expected_vectors
