import datetime as dt
import unittest
from unittest.mock import patch

import numpy as np

from python_raytrace.indices import resolve_space_weather_indices


class IndicesTests(unittest.TestCase):
    def test_resolve_indices_uses_pymsis_when_not_overridden(self) -> None:
        when = dt.datetime(2020, 1, 15, 12, 0, 0)
        with patch("python_raytrace.indices.pymsis.utils.get_f107_ap") as mock_get, patch(
            "python_raytrace.indices.refresh_pymsis_indices"
        ) as mock_refresh:
            mock_get.return_value = (
                np.array([120.0], dtype=float),
                np.array([110.0], dtype=float),
                np.array([[8.0, 9.0, 7.0, 6.0, 5.0, 4.5, 4.0]], dtype=float),
            )
            indices = resolve_space_weather_indices(when)

        self.assertEqual(indices.f107, 120.0)
        self.assertEqual(indices.f107a, 110.0)
        self.assertEqual(indices.ap_daily, 8.0)
        self.assertEqual(indices.source, "pymsis")
        self.assertEqual(indices.ap_vector.shape, (7,))
        mock_refresh.assert_not_called()

    def test_manual_overrides_replace_fetched_values(self) -> None:
        when = dt.datetime(2020, 1, 15, 12, 0, 0)
        with patch("python_raytrace.indices.pymsis.utils.get_f107_ap") as mock_get:
            mock_get.return_value = (
                np.array([120.0], dtype=float),
                np.array([110.0], dtype=float),
                np.array([[8.0, 9.0, 7.0, 6.0, 5.0, 4.5, 4.0]], dtype=float),
            )
            indices = resolve_space_weather_indices(when, f107=140.0, ap_daily=12.0)

        self.assertEqual(indices.f107, 140.0)
        self.assertEqual(indices.f107a, 140.0)
        self.assertEqual(indices.ap_daily, 12.0)
        self.assertEqual(indices.ap_vector[0], 12.0)
        self.assertEqual(indices.source, "manual")


if __name__ == "__main__":
    unittest.main()
