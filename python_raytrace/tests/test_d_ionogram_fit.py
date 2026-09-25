from __future__ import annotations

import unittest

import numpy as np

from reports.fit_d_ionogram import Ionogram, estimate_height_shift, score


FREQUENCIES = np.arange(2.0, 10.0001, .1)
SETTINGS = ("adaptive", "equal_area_guarded", .5, 4, 1000.0)


def ionogram(records: list[list[float]]) -> Ionogram:
    return Ionogram(FREQUENCIES, np.asarray(records, dtype=float), None, None, SETTINGS)


class DIonogramFitTests(unittest.TestCase):
    def test_all_accepted_paths_affect_score(self) -> None:
        observed = ionogram([[10, 1, 600], [10, 1, 710], [10, -1, 630]])
        incomplete = ionogram([[10, 1, 600], [10, -1, 630]])
        self.assertEqual(score(observed, observed)["total"], 0.0)
        self.assertGreater(score(observed, incomplete)["return_distance"], .1)

    def test_height_seed_uses_two_way_range(self) -> None:
        reference = ionogram([[i, mode, 700 + 3*i] for mode in (1, -1) for i in range(10, 30)])
        shifted = ionogram([[i, mode, 676 + 3*i] for mode in (1, -1) for i in range(10, 30)])
        self.assertAlmostEqual(estimate_height_shift(shifted, reference), 12.0)

    def test_nose_at_sweep_boundary_is_censored(self) -> None:
        a = ionogram([[80, 1, 600], [80, -1, 620]])
        b = ionogram([[80, 1, 610], [80, -1, 630]])
        self.assertEqual(score(a, b)["nose"], 0.0)


if __name__ == "__main__":
    unittest.main()
