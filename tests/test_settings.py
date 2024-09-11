from stjames import Constraint, Mode, OptimizationSettings, Settings


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
