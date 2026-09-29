from __future__ import annotations

import unittest

import numpy as np

from reports.ionogram_metrics import Ionogram, score


FREQUENCIES = np.arange(2.0, 10.0001, .1)
SETTINGS = ("adaptive", "equal_area_guarded", .5, 4, 1000.0)


def ionogram(records: list[list[float]]) -> Ionogram:
    return Ionogram(FREQUENCIES, np.asarray(records, dtype=float), None, None, SETTINGS)


class IonogramMetricsTests(unittest.TestCase):
    def test_all_accepted_paths_affect_score(self) -> None:
        observed = ionogram([[10, 1, 600], [10, 1, 710], [10, -1, 630]])
        incomplete = ionogram([[10, 1, 600], [10, -1, 630]])
        self.assertEqual(score(observed, observed)["total"], 0.0)
        self.assertGreater(score(observed, incomplete)["return_distance"], .1)

    def test_nose_at_sweep_boundary_is_censored(self) -> None:
        a = ionogram([[80, 1, 600], [80, -1, 620]])
        b = ionogram([[80, 1, 610], [80, -1, 630]])
        self.assertEqual(score(a, b)["nose"], 0.0)


if __name__ == "__main__":
    unittest.main()
