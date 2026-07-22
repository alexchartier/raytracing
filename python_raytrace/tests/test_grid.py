import datetime as dt
from pathlib import Path
import unittest

import numpy as np

from python_raytrace.geometry import GeoPoint
from python_raytrace.grid import IonosphereGrid, build_pyiri_grid
from python_raytrace.tracer import PointToPointRayTracer, RayTrace


class GridTests(unittest.TestCase):
    @unittest.skipUnless(
        Path("python_raytrace/_lib/libiri2020_bridge.dylib").exists() or Path("python_raytrace/_lib/libiri2020_bridge.so").exists(),
        "IRI2020 bridge library has not been built",
    )
    def test_pyiri_grid_builds_expected_shapes(self) -> None:
        grid = build_pyiri_grid(
            dt.datetime(2020, 1, 15, 12, 0, 0),
            GeoPoint(-20.0, 170.0, 0.0),
            GeoPoint(-18.0, -178.0, 0.0),
            f107=100.0,
            ap_daily=4.0,
            alt_min_km=70.0,
            alt_max_km=130.0,
            alt_step_km=10.0,
            lat_step_deg=2.0,
            lon_step_deg=2.0,
            lat_margin_deg=1.0,
            lon_margin_deg=1.0,
            d_region_model="fpt2018",
        )
        self.assertEqual(
            grid.iono_en_grid.shape,
            (len(grid.latitudes_deg), len(grid.longitudes_deg), len(grid.altitudes_km)),
        )
        self.assertEqual(grid.Bx.shape, grid.iono_en_grid.shape)
        self.assertEqual(grid.By.shape, grid.iono_en_grid.shape)
        self.assertEqual(grid.Bz.shape, grid.iono_en_grid.shape)
        self.assertGreater(grid.iono_en_grid.max(), 0.0)
        self.assertGreater(grid.collision_freq.max(), 0.0)
        self.assertGreaterEqual(grid.iono_grid_parms[3], -180.0)
        self.assertLessEqual(grid.iono_grid_parms[3], 180.0)
        self.assertIsNotNone(grid.electron_temp_k)
        self.assertIsNotNone(grid.ion_temp_k)
        self.assertIsNotNone(grid.neutral_temp_k)
        self.assertEqual(grid.collision_freq.shape, grid.iono_en_grid.shape)
        self.assertEqual(grid.metadata["source_grid"], "global_pyiri_then_regional_subset")
        self.assertGreater(grid.metadata["source_lon_count"], len(grid.longitudes_deg))

    def test_solver_homes_with_fake_backend(self) -> None:
        class FakeBackend:
            def trace(self, origin, elevations_deg, bearings_deg, freqs_mhz, ox_mode, nhops, tol, *, grid=None, state_vector=None):
                rays = []
                for elev, bear, freq in zip(elevations_deg, bearings_deg, freqs_mhz):
                    end_lat = 0.0 + 0.1 * (bear - 90.0)
                    end_lon = 1.0 + 0.1 * (elev - 40.0)
                    path = {
                        "initial_elev": float(elev),
                        "initial_bearing": float(bear),
                        "frequency": float(freq),
                        "lat": np.array([origin.lat_deg, end_lat], dtype=float),
                        "lon": np.array([origin.lon_deg, end_lon], dtype=float),
                        "height": np.array([origin.alt_km, 0.0], dtype=float),
                        "group_range": np.array([0.0, 100.0], dtype=float),
                        "geometric_distance": np.array([0.0, 100.0], dtype=float),
                        "absorption": np.array([0.0, 0.0], dtype=float),
                    }
                    rays.append(RayTrace(summary={}, path=path, state={}))
                return rays

        grid = IonosphereGrid(
            latitudes_deg=np.array([-1.0, 1.0], dtype=float),
            longitudes_deg=np.array([0.0, 2.0], dtype=float),
            altitudes_km=np.array([90.0, 100.0], dtype=float),
            iono_en_grid=np.zeros((2, 2, 2), dtype=float),
            iono_en_grid_5=np.zeros((2, 2, 2), dtype=float),
            collision_freq=np.zeros((2, 2, 2), dtype=float),
            iono_grid_parms=[-1.0, 2.0, 2.0, 0.0, 2.0, 2.0, 90.0, 10.0, 2.0],
            Bx=np.zeros((2, 2, 2), dtype=float),
            By=np.zeros((2, 2, 2), dtype=float),
            Bz=np.zeros((2, 2, 2), dtype=float),
            geomag_grid_parms=[-1.0, 2.0, 2.0, 0.0, 2.0, 2.0, 90.0, 10.0, 2.0],
        )
        tracer = PointToPointRayTracer(backend=FakeBackend())
        ray = tracer.trace_link(
            tx=GeoPoint(0.0, 0.0, 0.0),
            rx=GeoPoint(0.0, 1.0, 0.0),
            frequency_mhz=5.0,
            grid=grid,
            homing_tolerance_m=100.0,
            nhops=1,
        )
        self.assertTrue(ray.home)
        self.assertAlmostEqual(ray.launch_bearing_deg, 90.0, delta=1e-3)
        self.assertAlmostEqual(ray.launch_elevation_deg, 40.0, delta=1e-3)


if __name__ == "__main__":
    unittest.main()
