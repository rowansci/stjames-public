from pydantic import ValidationError
from pytest import raises

from stjames import Constraint
from stjames.conformers import ConformerGenMixin, ConformerSearchMixin, ETKDGSettings, iMTDSettings, iMTDSpeeds
from stjames.method import Method


def test_etkdg() -> None:
    settings = ETKDGSettings()

    assert settings.num_initial_confs == 300

    with raises(ValidationError, match="ETKDG does not support NCI"):
        ETKDGSettings(nci=True)

    with raises(ValidationError, match="ETKDG does not support constraints"):
        ETKDGSettings(constraints=[Constraint(constraint_type="bond", atoms=[1, 2])])


def test_imtdgc() -> None:
    settings = iMTDSettings()

    assert settings.speed == iMTDSpeeds.QUICK
    assert not settings.reopt
    assert settings.mtd_method == Method.GFN_FF


def test_conformer_gen_mixin() -> None:
    settings = ConformerGenMixin(conf_gen_settings=ETKDGSettings(num_initial_confs=150))
    assert settings.conf_gen_settings == ETKDGSettings(num_initial_confs=150)


def test_conformer_search_mixin() -> None:
    settings = ConformerSearchMixin(frequencies=True, solvent="water", conf_gen_settings=ETKDGSettings())

    assert settings.solvent == "water"
    assert settings.frequencies
