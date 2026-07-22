import datetime as dt
import tempfile
import unittest
from pathlib import Path

from python_raytrace.collm_skiymet_timeseries import (
    CollmSkiymetRadar,
    HookeWaveModel,
    TimeSeriesConfig,
    _resolve_wave_bearing_deg,
    build_cpa_swath_snapshot,
    build_fof2_map_snapshot,
    predict_collm_skiymet_time_series,
)


TLE_TEXT = """ISS (ZARYA)
1 25544U 98067A   26141.16510469  .00005835  00000+0  11282-3 0  9993
2 25544  51.6328  73.8715 0007528  81.3651 278.8190 15.49291753567564
"""


class CollmSkiymetTimeSeriesTests(unittest.TestCase):
    def test_wave_bearing_offset_applies_to_fallback_track_bearing(self) -> None:
        self.assertAlmostEqual(
            _resolve_wave_bearing_deg(HookeWaveModel(bearing_deg=None, bearing_offset_deg=-45.0), 120.0),
            75.0,
        )
        self.assertAlmostEqual(
            _resolve_wave_bearing_deg(HookeWaveModel(bearing_deg=30.0, bearing_offset_deg=15.0), 120.0),
            45.0,
        )

    def test_prediction_produces_folded_range_and_doppler(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tle_path = Path(tmpdir) / "iss.tle"
            tle_path.write_text(TLE_TEXT)
            series = predict_collm_skiymet_time_series(
                TimeSeriesConfig(
                    tle_file=tle_path,
                    start_utc=dt.datetime(2026, 6, 1, 15, 58, 35, tzinfo=dt.timezone.utc),
                    seconds=4.0,
                    step_seconds=2.0,
                    f107=100.0,
                    f107a=100.0,
                )
            )

        radar = CollmSkiymetRadar()
        self.assertEqual(len(series.points), 3)
        self.assertEqual(series.image_snr_db.shape[1], len(series.points))
        self.assertGreater(series.image_snr_db.shape[0], 100)
        for point in series.points:
            self.assertGreaterEqual(point.folded_range_km, 0.0)
            self.assertLessEqual(point.folded_range_km, radar.unambiguous_range_km)
            self.assertGreaterEqual(point.aliased_doppler_hz, -radar.doppler_nyquist_hz)
            self.assertLessEqual(point.aliased_doppler_hz, radar.doppler_nyquist_hz)
            self.assertTrue(point.home)
            self.assertIsNotNone(point.launch_bearing_deg)
            self.assertIsNotNone(point.launch_elevation_deg)

    def test_fof2_snapshot_contains_finite_grid_and_track(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tle_path = Path(tmpdir) / "iss.tle"
            tle_path.write_text(TLE_TEXT)
            series = predict_collm_skiymet_time_series(
                TimeSeriesConfig(
                    tle_file=tle_path,
                    start_utc=dt.datetime(2026, 6, 1, 15, 58, 35, tzinfo=dt.timezone.utc),
                    seconds=4.0,
                    step_seconds=2.0,
                    f107=100.0,
                    f107a=100.0,
                )
            )
            snapshot = build_fof2_map_snapshot(
                series,
                lat_step_deg=2.0,
                lon_step_deg=2.0,
                margin_deg=4.0,
                alt_min_km=180.0,
                alt_max_km=420.0,
                alt_step_km=20.0,
            )

        self.assertGreater(snapshot.latitudes_deg.size, 2)
        self.assertGreaterEqual(snapshot.longitudes_deg.size, 2)
        self.assertEqual(snapshot.fof2_mhz.shape, (snapshot.latitudes_deg.size, snapshot.longitudes_deg.size))
        self.assertTrue((snapshot.fof2_mhz > 0.0).any())
        self.assertEqual(snapshot.track_latitudes_deg.size, len(series.points))
        self.assertEqual(snapshot.track_longitudes_deg.size, len(series.points))

    def test_cpa_swath_contains_plasma_slice_and_ray(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tle_path = Path(tmpdir) / "iss.tle"
            tle_path.write_text(TLE_TEXT)
            series = predict_collm_skiymet_time_series(
                TimeSeriesConfig(
                    tle_file=tle_path,
                    start_utc=dt.datetime(2026, 6, 1, 15, 58, 35, tzinfo=dt.timezone.utc),
                    seconds=4.0,
                    step_seconds=2.0,
                    f107=100.0,
                    f107a=100.0,
                )
            )
            snapshot = build_cpa_swath_snapshot(
                series,
                waypoint_count=64,
                alt_min_km=0.0,
                alt_max_km=450.0,
                alt_step_km=10.0,
            )

        self.assertGreater(snapshot.along_track_km.size, 10)
        self.assertGreater(snapshot.altitudes_km.size, 10)
        self.assertEqual(snapshot.plasma_freq_mhz.shape, (snapshot.altitudes_km.size, snapshot.along_track_km.size))
        self.assertTrue((snapshot.plasma_freq_mhz > 0.0).any())
        self.assertGreater(snapshot.ray_along_track_km.size, 2)
        self.assertEqual(snapshot.ray_along_track_km.size, snapshot.ray_altitudes_km.size)
        self.assertGreaterEqual(float(snapshot.ray_along_track_km[-1]), float(snapshot.ray_along_track_km[0]))
        self.assertLessEqual(float(snapshot.ray_altitudes_km[0]), 5.0)

    def test_ephemeris_shift_and_hooke_wave_change_prediction(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            tle_path = Path(tmpdir) / "iss.tle"
            tle_path.write_text(TLE_TEXT)
            base = TimeSeriesConfig(
                tle_file=tle_path,
                start_utc=dt.datetime(2026, 6, 1, 15, 58, 35, tzinfo=dt.timezone.utc),
                seconds=4.0,
                step_seconds=2.0,
                f107=100.0,
                f107a=100.0,
                ephemeris_time_shift_seconds=0.0,
                hooke_wave=HookeWaveModel(amplitude_fraction=0.0, snr_coupling_db=0.0),
            )
            shifted = TimeSeriesConfig(
                tle_file=tle_path,
                start_utc=base.start_utc,
                seconds=base.seconds,
                step_seconds=base.step_seconds,
                f107=base.f107,
                f107a=base.f107a,
                ephemeris_time_shift_seconds=-3.0,
                hooke_wave=base.hooke_wave,
            )
            waved = TimeSeriesConfig(
                tle_file=tle_path,
                start_utc=base.start_utc,
                seconds=base.seconds,
                step_seconds=base.step_seconds,
                f107=base.f107,
                f107a=base.f107a,
                ephemeris_time_shift_seconds=base.ephemeris_time_shift_seconds,
                hooke_wave=HookeWaveModel(amplitude_fraction=0.16, horizontal_wavelength_km=100.0, period_seconds=20.0),
            )
            base_series = predict_collm_skiymet_time_series(base)
            shifted_series = predict_collm_skiymet_time_series(shifted)
            waved_series = predict_collm_skiymet_time_series(waved)

        base_ranges = [point.absolute_range_km for point in base_series.points]
        shifted_ranges = [point.absolute_range_km for point in shifted_series.points]
        base_snrs = [point.predicted_peak_snr_db for point in base_series.points]
        waved_snrs = [point.predicted_peak_snr_db for point in waved_series.points]
        self.assertNotEqual(base_ranges, shifted_ranges)
        self.assertTrue(any(abs(a - b) > 1e-6 for a, b in zip(base_snrs, waved_snrs)))


if __name__ == "__main__":
    unittest.main()
