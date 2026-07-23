from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable

import numpy as np


WGS84_A_M = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)


@dataclass(frozen=True)
class GeoPoint:
    lat_deg: float
    lon_deg: float
    alt_km: float = 0.0


@dataclass(frozen=True)
class RayDistance:
    distance_m: float
    closest_point_ecef_m: np.ndarray | None
    segment_index: int | None
    group_path_km: float | None
    geometric_path_km: float | None
    total_absorption_db: float | None


def wrap_longitude(lon_deg: float, minimum: float = -180.0) -> float:
    return ((lon_deg - minimum) % 360.0) + minimum


def wrap_longitudes(lons_deg: Iterable[float], minimum: float = -180.0) -> np.ndarray:
    lons = np.asarray(list(lons_deg), dtype=float)
    return ((lons - minimum) % 360.0) + minimum


def coerce_longitude_for_grid(lon_deg: float, grid_lons_deg: np.ndarray) -> float:
    candidates = np.array([lon_deg - 360.0, lon_deg, lon_deg + 360.0], dtype=float)
    center = 0.5 * (float(grid_lons_deg[0]) + float(grid_lons_deg[-1]))
    return float(candidates[np.argmin(np.abs(candidates - center))])


def llh_to_ecef(lat_deg: np.ndarray | float, lon_deg: np.ndarray | float, height_m: np.ndarray | float) -> np.ndarray:
    lat = np.deg2rad(np.asarray(lat_deg, dtype=float))
    lon = np.deg2rad(np.asarray(lon_deg, dtype=float))
    height = np.asarray(height_m, dtype=float)

    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)
    cos_lon = np.cos(lon)
    sin_lon = np.sin(lon)

    prime_vertical = WGS84_A_M / np.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
    x = (prime_vertical + height) * cos_lat * cos_lon
    y = (prime_vertical + height) * cos_lat * sin_lon
    z = (prime_vertical * (1.0 - WGS84_E2) + height) * sin_lat
    return np.stack((x, y, z), axis=-1)


def enu_to_ecef(east: np.ndarray | float, north: np.ndarray | float, up: np.ndarray | float,
                lat_deg: np.ndarray | float, lon_deg: np.ndarray | float) -> np.ndarray:
    east = np.asarray(east, dtype=float)
    north = np.asarray(north, dtype=float)
    up = np.asarray(up, dtype=float)
    lat = np.deg2rad(np.asarray(lat_deg, dtype=float))
    lon = np.deg2rad(np.asarray(lon_deg, dtype=float))

    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)
    sin_lon = np.sin(lon)
    cos_lon = np.cos(lon)

    x = -sin_lon * east - sin_lat * cos_lon * north + cos_lat * cos_lon * up
    y = cos_lon * east - sin_lat * sin_lon * north + cos_lat * sin_lon * up
    z = cos_lat * north + sin_lat * up
    return np.stack((x, y, z), axis=-1)


def relaz_to_ecef_unit(elevation_deg: float, azimuth_deg: float, lat_deg: float, lon_deg: float) -> np.ndarray:
    elev = math.radians(elevation_deg)
    az = math.radians(azimuth_deg)
    east = math.cos(elev) * math.sin(az)
    north = math.cos(elev) * math.cos(az)
    up = math.sin(elev)
    return enu_to_ecef(east, north, up, lat_deg, lon_deg)


def vector_angle_deg(vec_a: np.ndarray, vec_b: np.ndarray) -> float:
    norm_a = float(np.linalg.norm(vec_a))
    norm_b = float(np.linalg.norm(vec_b))
    if norm_a == 0.0 or norm_b == 0.0:
        return 0.0
    dot = float(np.dot(vec_a, vec_b) / (norm_a * norm_b))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot))))


def initial_bearing_deg(start: GeoPoint, end: GeoPoint) -> float:
    lat1 = math.radians(start.lat_deg)
    lat2 = math.radians(end.lat_deg)
    dlon = math.radians(wrap_longitude(end.lon_deg - start.lon_deg))
    y = math.sin(dlon) * math.cos(lat2)
    x = math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(dlon)
    return math.degrees(math.atan2(y, x))


def destination_point(start: GeoPoint, bearing_deg: float, distance_km: float, alt_km: float | None = None) -> GeoPoint:
    radius_km = 6371.0088
    angular_distance = float(distance_km) / radius_km
    lat1 = math.radians(start.lat_deg)
    lon1 = math.radians(start.lon_deg)
    bearing = math.radians(bearing_deg)

    lat2 = math.asin(
        math.sin(lat1) * math.cos(angular_distance)
        + math.cos(lat1) * math.sin(angular_distance) * math.cos(bearing)
    )
    lon2 = lon1 + math.atan2(
        math.sin(bearing) * math.sin(angular_distance) * math.cos(lat1),
        math.cos(angular_distance) - math.sin(lat1) * math.sin(lat2),
    )
    return GeoPoint(
        lat_deg=math.degrees(lat2),
        lon_deg=wrap_longitude(math.degrees(lon2)),
        alt_km=start.alt_km if alt_km is None else float(alt_km),
    )


def _unit_sphere_vector(lat_deg: float, lon_deg: float) -> np.ndarray:
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    return np.array([
        math.cos(lat) * math.cos(lon),
        math.cos(lat) * math.sin(lon),
        math.sin(lat),
    ], dtype=float)


def great_circle_waypoints(start: GeoPoint, end: GeoPoint, count: int = 128) -> tuple[np.ndarray, np.ndarray]:
    count = max(2, int(count))
    p0 = _unit_sphere_vector(start.lat_deg, start.lon_deg)
    p1 = _unit_sphere_vector(end.lat_deg, end.lon_deg)
    dot = max(-1.0, min(1.0, float(np.dot(p0, p1))))
    omega = math.acos(dot)
    t = np.linspace(0.0, 1.0, count)

    if omega < 1e-9:
        pts = np.outer(1.0 - t, p0) + np.outer(t, p1)
    else:
        sin_omega = math.sin(omega)
        pts = (
            np.sin((1.0 - t) * omega)[:, None] / sin_omega * p0[None, :]
            + np.sin(t * omega)[:, None] / sin_omega * p1[None, :]
        )

    pts /= np.linalg.norm(pts, axis=1)[:, None]
    lats = np.degrees(np.arctan2(pts[:, 2], np.hypot(pts[:, 0], pts[:, 1])))
    lons = np.degrees(np.arctan2(pts[:, 1], pts[:, 0]))
    return lats, lons


def make_regional_lat_lon_grids(start: GeoPoint, end: GeoPoint,
                                lat_step_deg: float, lon_step_deg: float,
                                lat_margin_deg: float = 5.0, lon_margin_deg: float = 5.0,
                                waypoint_count: int = 128) -> tuple[np.ndarray, np.ndarray]:
    path_lats, path_lons = great_circle_waypoints(start, end, count=waypoint_count)
    path_lons_unwrapped = np.degrees(np.unwrap(np.deg2rad(path_lons)))

    lat_min = max(-90.0, float(np.min(path_lats)) - lat_margin_deg)
    lat_max = min(90.0, float(np.max(path_lats)) + lat_margin_deg)
    lon_min = float(np.min(path_lons_unwrapped)) - lon_margin_deg
    lon_max = float(np.max(path_lons_unwrapped)) + lon_margin_deg

    lat_start = math.floor(lat_min / lat_step_deg) * lat_step_deg
    lat_stop = math.ceil(lat_max / lat_step_deg) * lat_step_deg
    lon_start = math.floor(lon_min / lon_step_deg) * lon_step_deg
    lon_stop = math.ceil(lon_max / lon_step_deg) * lon_step_deg

    while lon_start < -180.0:
        lon_start += 360.0
        lon_stop += 360.0
    while lon_start > 180.0:
        lon_start -= 360.0
        lon_stop -= 360.0

    lats = np.arange(lat_start, lat_stop + 0.5 * lat_step_deg, lat_step_deg, dtype=float)
    lats = np.clip(lats, -90.0, 90.0)
    lats = np.unique(lats)
    lons = np.arange(lon_start, lon_stop + 0.5 * lon_step_deg, lon_step_deg, dtype=float)
    return lats, lons


def expects_reflection(tx_alt_km: float, rx_alt_km: float) -> bool:
    if tx_alt_km < 100.0 and rx_alt_km < 100.0:
        return True
    if tx_alt_km > 400.0 and rx_alt_km > 400.0:
        return True
    return False


def _point_to_segment_distance(point: np.ndarray, seg_start: np.ndarray, seg_end: np.ndarray) -> tuple[float, float, np.ndarray]:
    delta = seg_end - seg_start
    denom = float(np.dot(delta, delta))
    if denom == 0.0:
        return float(np.linalg.norm(point - seg_start)), 0.0, seg_start
    t = float(np.dot(point - seg_start, delta) / denom)
    t = max(0.0, min(1.0, t))
    closest = seg_start + t * delta
    return float(np.linalg.norm(point - closest)), t, closest


def ray_point_distance(ray_path: dict, rx: GeoPoint, *, expect_reflection: bool | None = None) -> RayDistance:
    heights = np.asarray(ray_path.get("height", []), dtype=float)
    lats = np.asarray(ray_path.get("lat", []), dtype=float)
    lons = np.asarray(ray_path.get("lon", []), dtype=float)
    valid = np.isfinite(heights) & np.isfinite(lats) & np.isfinite(lons) & (heights < 1e40)

    if valid.sum() < 2:
        return RayDistance(math.inf, None, None, None, None, None)

    heights = heights[valid]
    lats = lats[valid]
    lons = lons[valid]

    start_index = 0
    if expect_reflection is None:
        expect_reflection = expects_reflection(float(heights[0]), rx.alt_km)
    if expect_reflection:
        apogee_index = int(np.argmax(heights))
        if apogee_index == len(heights) - 1:
            return RayDistance(math.inf, None, None, None, None, None)
        if float(np.min(heights[apogee_index:])) > rx.alt_km:
            return RayDistance(math.inf, None, None, None, None, None)
        start_index = apogee_index

    ray_xyz = llh_to_ecef(lats, lons, heights * 1000.0)
    point_xyz = llh_to_ecef(rx.lat_deg, rx.lon_deg, rx.alt_km * 1000.0)
    if point_xyz.ndim > 1:
        point_xyz = point_xyz[0]

    group_range = np.asarray(ray_path.get("group_range", np.full_like(heights, np.nan)), dtype=float)
    geom_distance = np.asarray(ray_path.get("geometric_distance", np.full_like(heights, np.nan)), dtype=float)
    absorption = np.asarray(ray_path.get("absorption", np.full_like(heights, np.nan)), dtype=float)

    seg_start = ray_xyz[start_index:-1]
    seg_end = ray_xyz[start_index + 1:]
    delta = seg_end - seg_start
    denom = np.sum(delta * delta, axis=1)
    rel = point_xyz[None, :] - seg_start
    t = np.zeros_like(denom)
    valid_denom = denom > 0.0
    t[valid_denom] = np.sum(rel[valid_denom] * delta[valid_denom], axis=1) / denom[valid_denom]
    t = np.clip(t, 0.0, 1.0)
    closest = seg_start + t[:, None] * delta
    distances = np.linalg.norm(point_xyz[None, :] - closest, axis=1)
    if distances.size == 0:
        return RayDistance(math.inf, None, None, None, None, None)
    best_offset = int(np.argmin(distances))
    best_distance = float(distances[best_offset])
    best_t = float(t[best_offset])
    best_segment = start_index + best_offset
    best_point = closest[best_offset]

    if best_segment is None or best_distance > 1e6:
        return RayDistance(math.inf, None, None, None, None, None)

    def _interp(values: np.ndarray) -> float | None:
        if values.shape[0] != heights.shape[0]:
            return None
        v0 = float(values[best_segment])
        v1 = float(values[best_segment + 1])
        if not math.isfinite(v0) or not math.isfinite(v1):
            return None
        return v0 + best_t * (v1 - v0)

    return RayDistance(
        distance_m=best_distance,
        closest_point_ecef_m=best_point,
        segment_index=best_segment,
        group_path_km=_interp(group_range),
        geometric_path_km=_interp(geom_distance),
        total_absorption_db=_interp(absorption),
    )
