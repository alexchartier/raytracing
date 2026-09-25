import math
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from python_raytrace.multisat_topside_inverse_demo import (
    DEFAULT_AMPERE_FILE,
    CaseObservables,
    HomedRayReturn,
    TopsideInverseConfig,
    _accumulate_doppler,
    _accumulate_image,
    _canonical_launch_angles,
    _deduplicate_homed_returns,
    _dense_plot_xlim_mhz,
    _masked_visible_field,
    _observed_support_extents,
    _plot_case_score,
    _profile_returns_in_window,
    _resample_vertical_display,
    build_inverse_problem,
    dataset_cost,
    make_frequency_axis_mhz,
    simulate_dataset,
    solve_inverse_problem,
)
from python_raytrace.grid import IonosphereGrid, load_ionosphere_grid_netcdf, save_ionosphere_grid_netcdf


class MultisatTopsideInverseDemoTests(unittest.TestCase):
    def test_return_intensity_does_not_depend_on_homing_miss(self) -> None:
        returns = (
            HomedRayReturn(2.0, 1, object(), 1.0, 100.0, 0.0, 10.0),
            HomedRayReturn(2.0, 1, object(), 999.0, 100.0, 0.0, 20.0),
        )
        image = np.zeros((1, 1), dtype=float)
        numerator = np.zeros_like(image)
        denominator = np.zeros_like(image)
        frequencies = np.array([2.0])
        range_edges = np.array([99.0, 101.0])
        _accumulate_image(image, returns=returns, frequencies_mhz=frequencies,
                          range_edges_km=range_edges)
        _accumulate_doppler(numerator, denominator, returns=returns,
                            frequencies_mhz=frequencies, range_edges_km=range_edges)
        self.assertEqual(float(image[0, 0]), 2.0)
        self.assertEqual(float(denominator[0, 0]), 2.0)
        self.assertEqual(float(numerator[0, 0]), 30.0)

    def test_equivalent_nadir_angles_are_one_return(self) -> None:
        elevation, bearing = _canonical_launch_angles(-92.0, 214.7)
        self.assertAlmostEqual(elevation, -88.0)
        self.assertAlmostEqual(bearing, 34.7)

        def ray_return(elev: float, bear: float, range_km: float) -> HomedRayReturn:
            ray = SimpleNamespace(path={"initial_elev": elev, "initial_bearing": bear})
            return HomedRayReturn(4.1, 1, ray, 100.0, range_km, 0.0, 0.0)

        returns = _deduplicate_homed_returns((
            ray_return(-92.0, 214.7, 958.05),
            ray_return(-88.0, 34.7, 957.39),
            ray_return(-87.9, 34.7, 957.50),
            ray_return(-82.0, 34.7, 957.50),
        ), 64)
        self.assertEqual(len(returns), 3)

    def test_frequency_axis_generation(self) -> None:
        axis = make_frequency_axis_mhz(5.0, 5.3, 100.0)
        np.testing.assert_allclose(axis, np.array([5.0, 5.1, 5.2, 5.3], dtype=float))

    def test_dense_plot_xlim_follows_refracted_nose(self) -> None:
        frequencies = np.array([3.0, 3.1, 3.2, 3.3, 3.4, 3.5], dtype=float)
        ranges = np.array([700.0, 750.0, 800.0, 1400.0, 1450.0, 1500.0], dtype=float)
        image = np.zeros((frequencies.size, ranges.size), dtype=float)
        image[:, 1] = 1.0
        image[0, 4] = 0.7
        image[1, 5] = 0.8
        xmin, xmax = _dense_plot_xlim_mhz(frequencies, ranges, image)
        self.assertAlmostEqual(xmin, 2.95)
        self.assertLess(xmax, 3.9)
        self.assertGreater(xmax, 3.25)

    def test_plot_helpers_prefer_supported_doppler(self) -> None:
        empty = CaseObservables(
            name="empty",
            kind="vertical",
            total_image=np.zeros((2, 3), dtype=float),
            o_image=np.zeros((2, 3), dtype=float),
            x_image=np.zeros((2, 3), dtype=float),
            ox_split_image=np.zeros((2, 3), dtype=float),
            doppler_image_hz=np.zeros((2, 3), dtype=float),
            doppler_weight=np.zeros((2, 3), dtype=float),
        )
        with_doppler = CaseObservables(
            name="with_doppler",
            kind="oblique",
            total_image=np.array([[0.0, 0.8, 0.0], [0.0, 1.0, 0.0]], dtype=float),
            o_image=np.zeros((2, 3), dtype=float),
            x_image=np.zeros((2, 3), dtype=float),
            ox_split_image=np.zeros((2, 3), dtype=float),
            doppler_image_hz=np.array([[0.0, -12.0, 0.0], [0.0, 8.0, 0.0]], dtype=float),
            doppler_weight=np.array([[0.0, 0.8, 0.0], [0.0, 1.0, 0.0]], dtype=float),
        )
        self.assertGreater(_plot_case_score(with_doppler), _plot_case_score(empty))
        masked, visible = _masked_visible_field(with_doppler.doppler_image_hz, with_doppler.doppler_weight)
        self.assertEqual(int(np.sum(~masked.mask)), 2)
        np.testing.assert_allclose(np.sort(visible), np.array([-12.0, 8.0], dtype=float))

    def test_observed_support_extents_follow_visible_nose_and_min_support(self) -> None:
        frequencies = np.array([1.0, 1.1, 1.2, 1.3], dtype=float)
        range_edges = np.array([700.0, 703.0, 706.0, 709.0, 712.0], dtype=float)
        image = np.zeros((frequencies.size, range_edges.size - 1), dtype=float)
        image[0, 0] = 0.2
        image[1, 1] = 0.6
        image[3, 3] = 1.0
        fmin, fmax, rmin, rmax = _observed_support_extents(image, frequencies, range_edges)
        self.assertAlmostEqual(fmin, 1.0)
        self.assertAlmostEqual(fmax, 1.3)
        self.assertAlmostEqual(rmin, 700.0)
        self.assertAlmostEqual(rmax, 712.0)

    def test_resample_vertical_display_changes_only_y_sampling(self) -> None:
        altitudes = np.array([0.0, 6.0, 12.0], dtype=float)
        field = np.array(
            [
                [0.0, 10.0],
                [6.0, 16.0],
                [12.0, 22.0],
            ],
            dtype=float,
        )
        display_alts, resampled = _resample_vertical_display(field, altitudes, step_km=3.0)
        np.testing.assert_allclose(display_alts, np.array([0.0, 3.0, 6.0, 9.0, 12.0], dtype=float))
        self.assertEqual(resampled.shape, (5, 2))
        np.testing.assert_allclose(resampled[:, 0], np.array([0.0, 3.0, 6.0, 9.0, 12.0], dtype=float))
        np.testing.assert_allclose(resampled[:, 1], np.array([10.0, 13.0, 16.0, 19.0, 22.0], dtype=float))

    def test_profile_returns_in_window_applies_window_and_miss_limit(self) -> None:
        returns = (
            HomedRayReturn(2.0, 1, object(), 10_000.0, 1200.0, 0.0, 0.0),
            HomedRayReturn(2.0, 1, object(), 25_000.0, 1250.0, 0.0, 0.0),
            HomedRayReturn(2.1, 1, object(), 15_000.0, 700.0, 0.0, 0.0),
            HomedRayReturn(2.1, 1, object(), 12_000.0, 1300.0, 0.0, 0.0),
            HomedRayReturn(5.5, 1, object(), 10_000.0, 1400.0, 0.0, 0.0),
        )
        selected = _profile_returns_in_window(
            returns,
            range_centers_km=np.array([700.0, 725.0, 1200.0, 1225.0, 1250.0, 1300.0, 1400.0], dtype=float),
            frequency_limits_mhz=(1.0, 5.0),
            max_miss_m=10_000.0,
        )
        self.assertEqual(len(selected), 1)
        self.assertAlmostEqual(selected[0].frequency_mhz, 2.0)
        self.assertAlmostEqual(selected[0].group_range_km, 1200.0)

    def test_grid_netcdf_roundtrip_supports_species_first_tensor(self) -> None:
        lat = np.array([0.0, 1.0], dtype=float)
        lon = np.array([10.0, 12.0], dtype=float)
        alt = np.array([90.0, 95.0, 100.0], dtype=float)
        shape = (lat.size, lon.size, alt.size)
        base = np.arange(np.prod(shape), dtype=float).reshape(shape)
        grid = IonosphereGrid(
            latitudes_deg=lat,
            longitudes_deg=lon,
            altitudes_km=alt,
            iono_en_grid=base,
            iono_en_grid_5=base + 1.0,
            collision_freq=base + 2.0,
            iono_grid_parms=[0.0] * 9,
            Bx=base + 3.0,
            By=base + 4.0,
            Bz=base + 5.0,
            geomag_grid_parms=[0.0] * 9,
            neutral_species_cm3=np.zeros((7,) + shape, dtype=float),
        )
        path = Path(tempfile.gettempdir()) / "grid_netcdf_roundtrip_test.nc"
        save_ionosphere_grid_netcdf(path, grid)
        loaded = load_ionosphere_grid_netcdf(path)
        self.assertEqual(loaded.neutral_species_cm3.shape, (7, lat.size, lon.size, alt.size))
        np.testing.assert_allclose(loaded.iono_en_grid, grid.iono_en_grid)
        from python_raytrace.grid import extract_ionosphere_subgrid

        subset = extract_ionosphere_subgrid(loaded, lat[1:], lon[:1], alt[1:])
        self.assertEqual(subset.neutral_species_cm3.shape, (7, 1, 1, 2))
        np.testing.assert_allclose(subset.neutral_species_cm3, grid.neutral_species_cm3[:, 1:, :1, 1:])

    @unittest.skipUnless(DEFAULT_AMPERE_FILE.exists(), "AMPERE test file is required for the topside inverse smoke test")
    def test_problem_and_solver_smoke(self) -> None:
        config = replace(
            TopsideInverseConfig(),
            frequencies_mhz=(5.0,),
            vertical_elevation_count=8,
            vertical_azimuth_step_deg=90.0,
            oblique_elevation_count=9,
            oblique_bearing_count=7,
            grid_lat_step_deg=2.0,
            grid_lon_step_deg=4.0,
            grid_alt_step_km=20.0,
            d_region_model="none",
            fit_parameter_names=("density_scale", "wave_phase_rad"),
            solver_maxiter=0,
            solver_popsize=4,
        )
        problem = build_inverse_problem(config)

        self.assertEqual(len(problem.cases), 8)
        observed = simulate_dataset(problem, config.truth_params)
        self.assertTrue(any(float(np.max(case.total_image)) > 0.0 for case in observed.cases))

        perturbed_params = replace(
            config.truth_params,
            density_scale=0.92,
            wave_phase_rad=config.truth_params.wave_phase_rad + 1.0,
        )
        truth_cost = dataset_cost(observed, observed)
        perturbed_cost = dataset_cost(observed, simulate_dataset(problem, perturbed_params))
        self.assertLess(truth_cost, perturbed_cost)

        fit = solve_inverse_problem(problem, observed)
        self.assertTrue(math.isfinite(fit.fitted_cost))
        self.assertGreaterEqual(fit.fitted_cost, fit.truth_cost)
        self.assertGreater(fit.evaluations, 0)


if __name__ == "__main__":
    unittest.main()
