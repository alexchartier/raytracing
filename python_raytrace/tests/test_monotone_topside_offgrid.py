from __future__ import annotations

import unittest

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from python_raytrace.general_field_inverse import LocalDensity
from python_raytrace.grid import IonosphereGrid
from python_raytrace.monotone_topside import monotone_topside_grid


class OffGridTrackTests(unittest.TestCase):
    def test_in_situ_anchor_between_grid_nodes(self) -> None:
        latitude = np.array([-2.0, 0.0, 2.0])
        longitude = np.array([0.0, 1.0])
        altitude = np.arange(100.0, 900.1, 20.0)
        peak = 200_000 * (1 + .02 * latitude[:, None] + .01 * longitude[None, :])
        distance = altitude - 300.0
        profile = peak[:, :, None] * np.exp(np.where(
            distance < 0, -.5 * (distance / 80.0) ** 2, -distance / 150.0))
        zero = np.zeros_like(profile)
        grid = IonosphereGrid(latitude, longitude, altitude, profile, profile.copy(),
                              zero, [0.0], zero, zero, zero, [0.0])
        track_lat = np.linspace(-.7, .7, 157)
        track_lon = .3
        observed = 20_000 + 1_500 * track_lat
        local = LocalDensity(np.column_stack((track_lat, np.full(len(track_lat), track_lon))),
                             800.0, observed, .005)
        result = monotone_topside_grid(grid, local, np.array([-.7, .7]),
                                       np.array([.72, .72]), np.array([.32, .32]),
                                       np.zeros(2))
        sampled = RegularGridInterpolator(
            (latitude, longitude, altitude), result.iono_en_grid)(
                np.column_stack((track_lat, np.full(len(track_lat), track_lon),
                                 np.full(len(track_lat), 800.0))))
        np.testing.assert_allclose(sampled, observed, rtol=.005)
        for row in result.iono_en_grid.reshape(-1, len(altitude)):
            self.assertTrue(np.all(np.diff(row[np.argmax(row):]) <= 1e-8 * row.max()))


if __name__ == "__main__":
    unittest.main()
