from pytest import raises

from stjames import Molecule
from stjames.workflows import MolecularDynamicsSettings, ThermodynamicEnsemble


def test_md_settings_raises(Mn: Molecule) -> None:
    MolecularDynamicsSettings()

    MolecularDynamicsSettings(temperature=100)

    MolecularDynamicsSettings(temperature=100, pressure=2)

    MolecularDynamicsSettings(temperature=None, ensemble=ThermodynamicEnsemble.NVE)

    with raises(ValueError):
        MolecularDynamicsSettings(ensemble=ThermodynamicEnsemble.NPT)

    with raises(ValueError):
        MolecularDynamicsSettings(temperature=None, ensemble=ThermodynamicEnsemble.NVT)
