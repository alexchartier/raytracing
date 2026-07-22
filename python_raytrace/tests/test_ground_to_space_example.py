import math
import unittest

from python_raytrace.ground_to_space_homing_example import (
    default_scenario,
    run_ground_to_space_demo,
)


class GroundToSpaceExampleTests(unittest.TestCase):
    def test_synthetic_demo_homes_onto_space_target(self) -> None:
        result = run_ground_to_space_demo("synthetic")
        scenario = default_scenario()

        self.assertTrue(result.ray.home)
        self.assertLess(result.ray.error_m, scenario.homing_tolerance_m)
        self.assertGreater(result.max_path_alt_km, scenario.rx.alt_km - 1.0)
        self.assertAlmostEqual(
            result.ray.launch_elevation_deg,
            result.line_of_sight_elevation_deg,
            delta=0.25,
        )
        self.assertAlmostEqual(
            result.ray.launch_bearing_deg,
            result.line_of_sight_bearing_deg,
            delta=0.25,
        )
        self.assertTrue(math.isfinite(result.slant_range_km))


if __name__ == "__main__":
    unittest.main()
