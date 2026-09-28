from types import SimpleNamespace

import numpy as np

from python_raytrace.geometry import GeoPoint
from python_raytrace.spacecraft_doppler import C_M_PER_S, spacecraft_doppler


def _ray(point: GeoPoint, tilt_deg: float) -> object:
    lat = np.array([point.lat_deg, point.lat_deg, point.lat_deg])
    lon = np.array([point.lon_deg, point.lon_deg, point.lon_deg])
    height = np.array([800.0, 600.0, 800.0])
    lat_rad = np.deg2rad(point.lat_deg)
    lon_rad = np.deg2rad(point.lon_deg)
    north = np.array([-np.sin(lat_rad) * np.cos(lon_rad),
                      -np.sin(lat_rad) * np.sin(lon_rad), np.cos(lat_rad)])
    up = np.array([np.cos(lat_rad) * np.cos(lon_rad),
                   np.cos(lat_rad) * np.sin(lon_rad), np.sin(lat_rad)])
    tilt = np.deg2rad(tilt_deg)
    departure = np.sin(tilt) * north - np.cos(tilt) * up
    arrival = -np.sin(tilt) * north + np.cos(tilt) * up
    momentum = np.stack((departure, np.zeros(3), arrival))
    return SimpleNamespace(
        path={"lat": lat, "lon": lon, "height": height},
        state={key: momentum[:, i] for i, key in
               enumerate(("dir_x", "dir_y", "dir_z"))},
    )


def test_vertical_return_has_zero_horizontal_spacecraft_doppler() -> None:
    point = GeoPoint(-50.0, 7.7, 800.0)
    result = spacecraft_doppler(_ray(point, 0.0), point, point, 5.0)
    assert abs(result.doppler_hz) < 1e-10


def test_projected_off_nadir_return_matches_two_way_doppler() -> None:
    point = GeoPoint(-50.0, 7.7, 800.0)
    result = spacecraft_doppler(_ray(point, 10.0), point, point, 5.0)
    expected = 2 * 5e6 * 8000 * np.sin(np.deg2rad(10.0)) / C_M_PER_S
    np.testing.assert_allclose(result.doppler_hz, expected, rtol=1e-12)
    np.testing.assert_allclose(result.launch_elevation_deg, -80.0, atol=1e-10)
    np.testing.assert_allclose(result.arrival_elevation_deg, 80.0, atol=1e-10)
