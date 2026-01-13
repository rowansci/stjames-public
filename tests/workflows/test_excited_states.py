"""Tests for excited state workflows."""

from stjames import Method, Settings, Task
from stjames.workflows import OmegaTuning, TDDFTSettings


def test_tddft_settings_basic() -> None:
    """Test basic TDDFTSettings construction."""
    settings = Settings(method=Method.B3LYP, basis_set="def2-SVP")
    tddft_settings = TDDFTSettings(settings=settings)

    assert tddft_settings.settings.method == Method.B3LYP
    assert tddft_settings.settings.basis_set
    assert tddft_settings.settings.basis_set.name == "def2-SVP"
    assert tddft_settings.tasks == {Task.ENERGY}
    assert tddft_settings.tda is True
    assert tddft_settings.num_excitations == 5
    assert tddft_settings.target_root is None
    assert tddft_settings.omega is None
    assert tddft_settings.settings_type == "TDDFTSettings"


def test_tddft_settings_custom() -> None:
    """Test TDDFTSettings with custom parameters."""
    settings = Settings(method=Method.WB97MD3BJ, basis_set="6-31G*")
    tddft_settings = TDDFTSettings(
        settings=settings,
        tasks={Task.ENERGY, Task.GRADIENT},
        tda=False,
        num_excitations=10,
        target_root=3,
        omega="koopmans",
    )

    assert tddft_settings.settings.method == Method.WB97MD3BJ
    assert tddft_settings.tasks == {Task.ENERGY, Task.GRADIENT}
    assert tddft_settings.tda is False
    assert tddft_settings.num_excitations == 10
    assert tddft_settings.target_root == 3
    assert tddft_settings.omega == OmegaTuning.KOOPMANS
