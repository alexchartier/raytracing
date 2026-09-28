"""Carrier Doppler from spacecraft motion through a frozen ionosphere.

PyLap's ray state ``dir_*`` is the phase refractive index times the local wave
normal in ECEF coordinates. Endpoint derivatives of phase path therefore give
the spacecraft contribution without retracing the same ray at displaced
positions. This differs from a derivative of group range.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .geometry import GeoPoint, llh_to_ecef

C_M_PER_S = 299_792_458.0


@dataclass(frozen=True)
class SpacecraftDoppler:
    doppler_hz: float
    launch_elevation_deg: float
    launch_bearing_deg: float
    arrival_elevation_deg: float
    arrival_bearing_deg: float
    receiver_miss_m: float


def _enu_basis(point: GeoPoint) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lat = math.radians(point.lat_deg)
    lon = math.radians(point.lon_deg)
    east = np.array([-math.sin(lon), math.cos(lon), 0.0])
    north = np.array([-math.sin(lat) * math.cos(lon),
                      -math.sin(lat) * math.sin(lon), math.cos(lat)])
    up = np.array([math.cos(lat) * math.cos(lon),
                   math.cos(lat) * math.sin(lon), math.sin(lat)])
    return east, north, up


def _angles(momentum: np.ndarray, point: GeoPoint) -> tuple[float, float]:
    east, north, up = _enu_basis(point)
    e, n, u = (float(np.dot(momentum, axis)) for axis in (east, north, up))
    return (math.degrees(math.atan2(u, math.hypot(e, n))),
            math.degrees(math.atan2(e, n)) % 360.0)


def spacecraft_doppler(ray: object, tx: GeoPoint, rx: GeoPoint,
                       frequency_mhz: float, speed_mps: float = 8000.0,
                       track_bearing_deg: float = 0.0) -> SpacecraftDoppler:
    """Predict Doppler for an Earth-fixed static medium and moving endpoints.

    ``track_bearing_deg`` is clockwise from north. For the monostatic pass the
    transmitter and receiver share this velocity; the echo's light-time motion
    is negligible at the current 1 km homing gate. The returned arrival angles
    describe propagation toward the receiver, not antenna pointing.
    """
    if speed_mps < 0 or not math.isfinite(speed_mps):
        raise ValueError("speed_mps must be finite and nonnegative")
    path = ray.path
    state = ray.state
    heights = np.asarray(path["height"], dtype=float)
    lats = np.asarray(path["lat"], dtype=float)
    lons = np.asarray(path["lon"], dtype=float)
    valid = (np.isfinite(heights) & np.isfinite(lats) & np.isfinite(lons)
             & (heights < 1e40))
    if np.count_nonzero(valid) < 3:
        raise ValueError("Ray path has too few valid points")
    points = llh_to_ecef(lats[valid], lons[valid], heights[valid] * 1000.0)
    momentum = np.column_stack([np.asarray(state[key], dtype=float)[valid]
                                for key in ("dir_x", "dir_y", "dir_z")])
    receiver = np.asarray(llh_to_ecef(rx.lat_deg, rx.lon_deg, rx.alt_km * 1000.0),
                          dtype=float).reshape(3)
    start = 1
    if tx.alt_km > 400 and rx.alt_km > 400:
        below = np.flatnonzero(heights[valid] < min(tx.alt_km, rx.alt_km) - 20.0)
        if len(below):
            start = int(below[0])
    if start >= len(points) - 1:
        raise ValueError("Ray does not have a return leg")
    segments = points[start + 1:] - points[start:-1]
    offsets = receiver - points[start:-1]
    lengths_squared = np.einsum("ij,ij->i", segments, segments)
    fractions = np.clip(np.divide(
        np.einsum("ij,ij->i", offsets, segments), lengths_squared,
        out=np.zeros_like(lengths_squared), where=lengths_squared > 0), 0, 1)
    closest = points[start:-1] + fractions[:, None] * segments
    distances = np.linalg.norm(closest - receiver, axis=1)
    local_index = int(np.argmin(distances))
    index = start + local_index
    fraction = float(fractions[local_index])
    departure = momentum[0]
    arrival = (1.0 - fraction) * momentum[index] + fraction * momentum[index + 1]
    east, north, _ = _enu_basis(tx)
    bearing = math.radians(track_bearing_deg)
    velocity = speed_mps * (math.sin(bearing) * east + math.cos(bearing) * north)
    # The two endpoint derivatives of optical phase path are -p_launch and
    # +p_arrival. Doppler is minus the time derivative of that phase path.
    doppler = frequency_mhz * 1e6 / C_M_PER_S * float(
        np.dot(velocity, departure - arrival))
    launch_elev, launch_bear = _angles(departure, tx)
    arrival_elev, arrival_bear = _angles(arrival, rx)
    return SpacecraftDoppler(
        doppler_hz=doppler,
        launch_elevation_deg=launch_elev,
        launch_bearing_deg=launch_bear,
        arrival_elevation_deg=arrival_elev,
        arrival_bearing_deg=arrival_bear,
        receiver_miss_m=float(distances[local_index]),
    )
