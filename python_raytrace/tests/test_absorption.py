import unittest

import numpy as np

from python_raytrace.absorption import effective_collision_frequency
from python_raytrace.iri2020_model import blend_low_altitude_density


class AbsorptionTests(unittest.TestCase):
    def test_effective_collision_frequency_is_positive(self) -> None:
        te = np.array([300.0, 500.0, 900.0], dtype=float)
        ti = np.array([280.0, 450.0, 850.0], dtype=float)
        ne = np.array([1e9, 5e10, 3e11], dtype=float)
        neutrals = np.zeros((7, 3), dtype=float)
        neutrals[2] = np.array([1e12, 5e11, 1e11], dtype=float)  # N2
        neutrals[3] = np.array([3e11, 2e11, 5e10], dtype=float)  # O2
        neutrals[1] = np.array([5e10, 4e10, 3e10], dtype=float)  # O

        coll = effective_collision_frequency(te, ti, ne, neutrals)
        self.assertEqual(coll.shape, ne.shape)
        self.assertTrue(np.all(coll > 0.0))

    def test_d_region_blend_prefers_low_altitude_profile(self) -> None:
        alts = np.array([80.0, 100.0, 120.0, 130.0, 140.0, 160.0], dtype=float)
        upper = np.array([1e7, 3e8, 2e9, 5e9, 1e10, 2e10], dtype=float)
        lower = np.array([1e9, 8e8, 5e8, 2e9, 8e9, 1.5e10], dtype=float)

        merged = blend_low_altitude_density(alts, upper, lower, blend_bottom_km=120.0, blend_top_km=140.0)
        self.assertAlmostEqual(merged[0], lower[0])
        self.assertAlmostEqual(merged[2], lower[2])
        self.assertAlmostEqual(merged[4], upper[4])
        self.assertGreater(merged[3], lower[3])
        self.assertLess(merged[3], upper[3])


if __name__ == "__main__":
    unittest.main()
