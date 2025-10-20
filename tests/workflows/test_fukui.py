from pytest import fixture

from stjames import Method, Molecule, Settings
from stjames.workflows import FukuiIndexWorkflow


@fixture
def water() -> Molecule:
    return Molecule.from_xyz("H 0 0 0\nO 0 0 1\nH 0 1 1")


def test_fukui_basic(water: Molecule) -> None:
    fi = FukuiIndexWorkflow(initial_molecule=water)
    assert fi.fukui_settings.method == Method.GFN1_XTB
    assert fi.opt_settings is None

    gfn2_settings = Settings(method=Method.GFN2_XTB)
    fj = FukuiIndexWorkflow(initial_molecule=water, opt_settings=gfn2_settings)
    assert fj.fukui_settings.method == Method.GFN1_XTB
    assert fj.opt_settings == gfn2_settings
