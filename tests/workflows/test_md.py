from pytest import raises

from stjames.workflows import MolecularDynamicsSettings, ThermodynamicEnsemble


def test_md_settings_raises() -> None:
    MolecularDynamicsSettings()

    MolecularDynamicsSettings(temperature=100)

    MolecularDynamicsSettings(temperature=100, pressure=2)

    with raises(ValueError):
        MolecularDynamicsSettings(ensemble=ThermodynamicEnsemble.NPT)
