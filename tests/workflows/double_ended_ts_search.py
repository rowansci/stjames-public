from pytest import raises

from stjames import Method, Settings
from stjames.molecule import Molecule
from stjames.optimization.freezing_string_method import FSMInterpolation, FSMOptimizationCoordinates, FSMSettings
from stjames.workflows.double_ended_ts_search import DoubleEndedTSSearchWorkflow


def test_setup() -> None:
    ammonia = Molecule.from_xyz("""\
N  0.0  0.0  0.3
H -0.5 -0.8 -0.1
H -0.4  0.8 -0.1
H  0.9  0.0 -0.1""")
    ammonia_flipped = Molecule.from_xyz("""\
N  0.0  0.0 -0.3
H -0.5 -0.8  0.1
H -0.4  0.8  0.1
H  0.9  0.0  0.1""")
    hcn = Molecule.from_xyz("""\
H 0 0 -1.1
C 0 0 0
N 0 0 1.2""")

    with raises(ValueError, match="Reactants and products must have the same number of atoms."):
        DoubleEndedTSSearchWorkflow(
            reactant=ammonia,
            product=hcn,
            calculation_settings=Settings(method=Method.GFN0_XTB),
            search_settings=FSMSettings(),
            optimize_inputs=True,
            optimize_ts=False,
        )

    DoubleEndedTSSearchWorkflow(
        reactant=ammonia,
        product=ammonia_flipped,
        calculation_settings=Settings(method=Method.HF3C),
        search_settings=FSMSettings(
            optimization_coordinates=FSMOptimizationCoordinates.REDUNDANT_INTERNAL_COORDINATES,
            interpolation_method=FSMInterpolation.LINEAR_SYNCHRONOUS_TRANSIT,
        ),
        optimize_inputs=False,
        optimize_ts=True,
    )
