import numpy as np

from reports.generate_lat_wave_ionogram import make_fan


def test_option_d_stays_default_and_option_a_uses_full_az_el_fan() -> None:
    default_config, default_elevations, default_bearings = make_fan()
    d_config, d_elevations, d_bearings = make_fan("D")
    a_config, a_elevations, a_bearings = make_fan("A")
    assert default_config == d_config
    np.testing.assert_array_equal(default_elevations, d_elevations)
    np.testing.assert_array_equal(default_bearings, d_bearings)
    assert d_config.vertical_fan_layout == "equal_area_guarded"
    assert d_config.vertical_guard_seed_limit == 4
    assert len(d_elevations) == 162
    assert a_config.vertical_fan_layout == "az_el"
    assert a_config.vertical_guard_seed_limit is None
    assert len(a_elevations) == len(a_bearings) == 288
