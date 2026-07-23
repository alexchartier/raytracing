import math
import unittest

from python_raytrace.ground_to_ground_homing_example import (
    build_electron_density_swath,
    default_scenario,
    run_ground_to_ground_demo,
)


class GroundToGroundExampleTests(unittest.TestCase):
    def test_default_demo_homes_onto_500_km_ground_target_in_both_modes(self) -> None:
        result = run_ground_to_ground_demo()
        scenario = default_scenario()

        self.assertTrue(result.ordinary_ray.home)
        self.assertTrue(result.extraordinary_ray.home)
        self.assertLess(result.ordinary_ray.error_m, scenario.homing_tolerance_m)
        self.assertLess(result.extraordinary_ray.error_m, scenario.homing_tolerance_m)
        self.assertAlmostEqual(result.great_circle_distance_km, 500.0, delta=0.5)
        self.assertAlmostEqual(result.great_circle_bearing_deg, scenario.bearing_deg, delta=0.5)
        self.assertGreater(result.ordinary_max_path_alt_km, scenario.tx.alt_km)
        self.assertGreater(result.extraordinary_max_path_alt_km, scenario.tx.alt_km)
        self.assertTrue(math.isfinite(result.ordinary_ray.launch_elevation_deg))
        self.assertTrue(math.isfinite(result.extraordinary_ray.launch_elevation_deg))

    def test_swath_contains_finite_density_and_ray_path(self) -> None:
        result = run_ground_to_ground_demo()
        swath = build_electron_density_swath(result, waypoint_count=64)

        self.assertEqual(swath.electron_density_m3.shape, (swath.altitudes_km.size, swath.along_track_km.size))
        self.assertTrue((swath.electron_density_m3 > 0.0).any())
        self.assertGreater(swath.ordinary_ray_along_track_km.size, 2)
        self.assertGreater(swath.extraordinary_ray_along_track_km.size, 2)
        self.assertEqual(swath.ordinary_ray_along_track_km.size, swath.ordinary_ray_altitudes_km.size)
        self.assertEqual(swath.extraordinary_ray_along_track_km.size, swath.extraordinary_ray_altitudes_km.size)


if __name__ == "__main__":
    unittest.main()
