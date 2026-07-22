import math
import unittest

import numpy as np

from python_raytrace.geometry import (
    GeoPoint,
    coerce_longitude_for_grid,
    make_regional_lat_lon_grids,
    ray_point_distance,
)


class GeometryTests(unittest.TestCase):
    def test_regional_grid_handles_dateline_crossing(self) -> None:
        tx = GeoPoint(10.0, 179.0, 0.0)
        rx = GeoPoint(12.0, -179.0, 0.0)
        lats, lons = make_regional_lat_lon_grids(
            tx,
            rx,
            lat_step_deg=1.0,
            lon_step_deg=1.0,
            lat_margin_deg=2.0,
            lon_margin_deg=2.0,
        )
        self.assertGreaterEqual(lons[0], -180.0)
        self.assertLessEqual(lons[0], 180.0)
        aligned = coerce_longitude_for_grid(rx.lon_deg, lons)
        self.assertGreaterEqual(aligned, lons[0] - 1.0)
        self.assertLessEqual(aligned, lons[-1] + 1.0)
        self.assertGreater(len(lats), 1)
        self.assertGreater(len(lons), 1)

    def test_ray_point_distance_interpolates_along_segment(self) -> None:
        ray_path = {
            "lat": np.array([0.0, 0.0], dtype=float),
            "lon": np.array([0.0, 1.0], dtype=float),
            "height": np.array([0.0, 0.0], dtype=float),
            "group_range": np.array([0.0, 100.0], dtype=float),
            "geometric_distance": np.array([0.0, 100.0], dtype=float),
            "absorption": np.array([0.0, 1.0], dtype=float),
        }
        rx = GeoPoint(0.0, 0.5, 0.0)
        distance = ray_point_distance(ray_path, rx)
        self.assertLess(distance.distance_m, 1e3)
        self.assertTrue(math.isfinite(distance.group_path_km))
        self.assertAlmostEqual(distance.group_path_km, 50.0, delta=5.0)


if __name__ == "__main__":
    unittest.main()
