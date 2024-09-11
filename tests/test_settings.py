from stjames import Constraint, Mode, OptimizationSettings, Settings


def test_set_mode_auto() -> None:
    Settings()
    assert Settings().mode == Mode.RAPID


def test_constraints_opt_settings() -> None:
    settings_rapid = Settings(mode=Mode.RAPID)
    settings_meticulous = Settings(mode=Mode.METICULOUS)

    cons = [Constraint(atoms=[1, 2], constraint_type="bond")]
    settings_careful_con = Settings(mode=Mode.CAREFUL, opt_settings=OptimizationSettings(constraints=cons))

    rap_opt_set = settings_rapid.opt_settings
    car_con_opt_set = settings_careful_con.opt_settings
    met_opt_set = settings_meticulous.opt_settings

    assert not rap_opt_set.constraints
    assert not met_opt_set.constraints
    assert car_con_opt_set.constraints == cons

    assert car_con_opt_set.energy_threshold == 5e-5
    # CAREFUL + constraints optimization settings == RAPID optimization settings
    assert car_con_opt_set.energy_threshold == rap_opt_set.energy_threshold
    assert car_con_opt_set.max_gradient_threshold == rap_opt_set.max_gradient_threshold
    assert car_con_opt_set.rms_gradient_threshold == rap_opt_set.rms_gradient_threshold
