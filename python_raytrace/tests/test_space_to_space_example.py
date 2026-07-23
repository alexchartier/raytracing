import math
import unittest
from dataclasses import replace

from python_raytrace.space_to_space_state_vector_example import (
    default_scenario,
    run_space_to_space_demo,
)


class SpaceToSpaceExampleTests(unittest.TestCase):
    def test_default_demo_homes_with_explicit_state_vectors(self) -> None:
        scenario = replace(
            default_scenario(),
            frequency_start_mhz=9.0,
            frequency_stop_mhz=10.0,
            frequency_step_khz=1000.0,
        )
        result = run_space_to_space_demo(scenario)

        self.assertEqual(len(result.sweep), 2)
        self.assertEqual([entry.frequency_mhz for entry in result.sweep], [9.0, 10.0])
        self.assertTrue(all(entry.ordinary.ray.home for entry in result.sweep))
        self.assertTrue(all(entry.extraordinary.ray.home for entry in result.sweep))
        self.assertTrue(all(entry.ordinary.ray.error_m < scenario.homing_tolerance_m for entry in result.sweep))
        self.assertTrue(all(entry.extraordinary.ray.error_m < scenario.homing_tolerance_m for entry in result.sweep))
        self.assertTrue(math.isfinite(result.line_of_sight_elevation_deg))
        self.assertLess(result.line_of_sight_elevation_deg, 0.0)
        self.assertTrue(all(entry.ordinary.seed_error_m > entry.ordinary.ray.error_m for entry in result.sweep))
        self.assertTrue(all(entry.extraordinary.seed_error_m > entry.extraordinary.ray.error_m for entry in result.sweep))


if __name__ == "__main__":
    unittest.main()
