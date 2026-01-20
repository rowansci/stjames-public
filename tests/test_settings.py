from pytest import raises

from stjames import Constraint, Method, Mode, OmegaTuning, OptimizationSettings, Settings, TDDFTSettings


def test_set_mode_auto() -> None:
    Settings()
    assert Settings().mode == Mode.RAPID


def test_opt_settings() -> None:
    settings_rapid = Settings(mode=Mode.RAPID)
    settings_meticulous = Settings(mode=Mode.METICULOUS)

    cons = [Constraint(atoms=[1, 2], constraint_type="bond")]
    settings_careful = Settings(mode=Mode.CAREFUL, opt_settings=OptimizationSettings(constraints=cons))

    rap_opt_set = settings_rapid.opt_settings
    car_opt_set = settings_careful.opt_settings
    met_opt_set = settings_meticulous.opt_settings

    assert not rap_opt_set.constraints
    assert not met_opt_set.constraints
    assert car_opt_set.constraints == cons

    assert rap_opt_set.energy_threshold == 5e-5
    assert rap_opt_set.max_gradient_threshold == 5e-3
    assert rap_opt_set.rms_gradient_threshold == 3.5e-3

    assert car_opt_set.energy_threshold == 1e-6
    assert car_opt_set.max_gradient_threshold == 9e-4
    assert car_opt_set.rms_gradient_threshold == 6e-4

    assert met_opt_set.energy_threshold == 1e-6
    assert met_opt_set.max_gradient_threshold == 3e-5
    assert met_opt_set.rms_gradient_threshold == 2e-5


def test_omega() -> None:
    with raises(ValueError, match="Omega tuning may only be specified for range-separated DFT functionals"):
        Settings(method=Method.B3LYP, omega="koopmans")

    Settings(method=Method.CAMB3LYP, omega="koopmans")
    Settings(method=Method.CAMB3LYP, omega=0.3)


def test_tddft_settings_basic() -> None:
    """Test basic TDDFTSettings construction."""
    settings = Settings(
        method=Method.WB97MD3BJ,
        basis_set="def2-SVP",
        omega="koopmans",
        excited_state_settings=TDDFTSettings(
            tda=False,
            num_excitations=8,
            target_root=2,
        ),
        engine="pyscf",
    )

    assert settings.method == Method.WB97MD3BJ
    assert settings.omega == OmegaTuning.KOOPMANS
    assert isinstance(settings.excited_state_settings, TDDFTSettings)
    assert settings.excited_state_settings.tda is False
    assert settings.excited_state_settings.num_excitations == 8
    assert settings.excited_state_settings.target_root == 2


def test_tddft_settings_custom() -> None:
    """Test TDDFTSettings with custom parameters."""
    settings = Settings(
        method=Method.B3LYP,
        basis_set="sto-3g",
        excited_state_settings=TDDFTSettings(
            tda=True,
            num_excitations=10,
            target_root=5,
        ),
        engine="gpu4pyscf",
    )

    assert settings.method == Method.B3LYP
    assert settings.omega is None
    assert isinstance(settings.excited_state_settings, TDDFTSettings)
    assert settings.excited_state_settings.tda is True
    assert settings.excited_state_settings.num_excitations == 10
    assert settings.excited_state_settings.target_root == 5


def test_settings_roundtrip() -> None:
    settings = Settings(
        method=Method.B3LYP,
        basis_set="sto-3g",
        excited_state_settings=TDDFTSettings(num_excitations=5),
        engine="pyscf",
    )

    data = settings.model_dump(mode="json")
    assert isinstance(data, dict)
    assert isinstance(data["excited_state_settings"], dict)

    stjames_settings = Settings.model_validate(data)

    assert stjames_settings.method == settings.method
    assert stjames_settings.basis_set == settings.basis_set
    assert stjames_settings.engine == settings.engine
    assert isinstance(stjames_settings.excited_state_settings, TDDFTSettings)
    assert stjames_settings.excited_state_settings == settings.excited_state_settings
