import numpy as np

from reports.oblique_wave_pass import FREQUENCIES, RUN, endpoints, separation_km


def test_saved_oblique_truth_contains_only_homed_reflections() -> None:
    for index in range(1, 21):
        tx, rx = endpoints(index)
        assert rx.lat_deg > tx.lat_deg
        np.testing.assert_allclose(separation_km(tx, rx), 600.0, atol=1e-6)
        path = RUN / "truth" / f"ionogram_{index:02d}.npz"
        with np.load(path, allow_pickle=False) as saved:
            records = np.asarray(saved["records"], dtype=float)
            assert str(saved["doppler_model"]) == "local transmitter and receiver tangent velocities"
            assert str(saved["method"]) == "oblique_adaptive"
            np.testing.assert_allclose(saved["frequencies_mhz"], FREQUENCIES)
            assert int(np.sum(saved["count_array"])) == len(records)
            assert len(saved["spacecraft_doppler_hz"]) == len(records)
            assert np.all(records[:, 3] <= 1000.0)
            assert np.all(records[:, 2] >= 700.0)
            assert np.all(saved["ray_minimum_altitude_km"] < 700.0)
            assert set(np.unique(records[:, 1])) == {-1.0, 1.0}
