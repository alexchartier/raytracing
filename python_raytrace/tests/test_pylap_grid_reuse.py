import unittest
from types import SimpleNamespace

import numpy as np

from python_raytrace.geometry import GeoPoint
from python_raytrace.tracer import PyLapRaytraceBackend


class PyLapGridReuseTests(unittest.TestCase):
    def setUp(self) -> None:
        PyLapRaytraceBackend._native_grid = None
        PyLapRaytraceBackend._native_func = None
        PyLapRaytraceBackend._reuse_unsupported_func = None

    tearDown = setUp

    def _grid(self):
        values = np.zeros((1, 1, 1), dtype=float)
        return SimpleNamespace(
            iono_en_grid=values, iono_en_grid_5=values,
            collision_freq=values, iono_grid_parms=[0.0] * 9,
            Bx=values, By=values, Bz=values, geomag_grid_parms=[0.0] * 9,
        )

    def _trace(self, backend, grid):
        return backend.trace(
            GeoPoint(0.0, 0.0, 100.0), [0.0], [0.0], [4.0], 1, 1,
            (1e-8, 0.005, 5.0), grid=grid,
            state_vector={"pos_x": np.array([1.0])},
        )

    def test_reuses_only_the_active_grid(self) -> None:
        lengths = []

        def raytrace(*args):
            lengths.append(len(args))
            return ([{}], [{}], [{}])

        backend = PyLapRaytraceBackend(raytrace, cache_native_grid=True)
        first_grid, second_grid = self._grid(), self._grid()
        self._trace(backend, first_grid)
        self._trace(backend, first_grid)
        self._trace(backend, second_grid)
        self._trace(backend, second_grid)
        self.assertEqual(lengths, [18, 10, 18, 10])

    def test_old_extension_falls_back_to_full_grid(self) -> None:
        lengths = []

        def raytrace(*args):
            lengths.append(len(args))
            if len(args) == 10:
                raise TypeError("argument 10 must be numpy.ndarray, not dict")
            return ([{}], [{}], [{}])

        backend = PyLapRaytraceBackend(raytrace, cache_native_grid=True)
        grid = self._grid()
        self._trace(backend, grid)
        self._trace(backend, grid)
        self._trace(backend, grid)
        self.assertEqual(lengths, [18, 10, 18, 18])

    def test_failed_grid_load_invalidates_previous_grid(self) -> None:
        lengths = []
        first_grid, second_grid = self._grid(), self._grid()

        def raytrace(*args):
            lengths.append(len(args))
            if len(args) == 18 and args[9] is second_grid.iono_en_grid:
                raise ValueError("grid load failed")
            return ([{}], [{}], [{}])

        backend = PyLapRaytraceBackend(raytrace, cache_native_grid=True)
        self._trace(backend, first_grid)
        with self.assertRaisesRegex(ValueError, "grid load failed"):
            self._trace(backend, second_grid)
        self._trace(backend, first_grid)
        self.assertEqual(lengths, [18, 18, 18])


if __name__ == "__main__":
    unittest.main()
