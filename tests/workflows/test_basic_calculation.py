from stjames.engine import Engine
from stjames.molecule import Molecule
from stjames.workflows import BasicCalculationWorkflow


def test_basic_calculation() -> None:
    """Test the basic calculation workflow."""
    mol = Molecule.from_xyz("He 0 0 0")
    workflow = BasicCalculationWorkflow(
        initial_molecule=mol,
        settings={"method": "GFN2-xTB"},
    )

    assert workflow.engine == Engine.XTB
    assert workflow.calculation_uuid is None
