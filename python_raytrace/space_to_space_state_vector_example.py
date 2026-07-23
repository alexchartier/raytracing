from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import griddata

if __package__:
    from .geometry import GeoPoint, coerce_longitude_for_grid, destination_point, great_circle_waypoints, llh_to_ecef, ray_point_distance, wrap_longitude
    from .grid import IonosphereGrid, build_pyiri_grid, load_ionosphere_grid, save_ionosphere_grid
    from .tracer import PointToPointRayTracer, PyLapImportError, RayTrace
else:
    from geometry import GeoPoint, coerce_longitude_for_grid, destination_point, great_circle_waypoints, llh_to_ecef, ray_point_distance, wrap_longitude
    from grid import IonosphereGrid, build_pyiri_grid, load_ionosphere_grid, save_ionosphere_grid
    from tracer import PointToPointRayTracer, PyLapImportError, RayTrace


@dataclass(frozen=True)
class SpaceToSpaceScenario:
    when: dt.datetime
    tx: GeoPoint
    distance_km: float = 500.0
    bearing_deg: float = 90.0
    rx_alt_km: float = 500.0
    frequency_start_mhz: float = 5.0
    frequency_stop_mhz: float = 15.0
    frequency_step_khz: float = 200.0
    f107: float = 120.0
    ap_daily: float = 8.0
    nhops: int = 1
    homing_tolerance_m: float = 50.0

    @property
    def frequencies_mhz(self) -> np.ndarray:
        step_mhz = max(self.frequency_step_khz / 1000.0, 1e-6)
        count = int(math.floor((self.frequency_stop_mhz - self.frequency_start_mhz) / step_mhz + 1e-9)) + 1
        values = self.frequency_start_mhz + step_mhz * np.arange(max(count, 1), dtype=float)
        return values[values <= self.frequency_stop_mhz + 1e-9]


@dataclass(frozen=True)
class ModeSolution:
    label: str
    ox_mode: int
    frequency_mhz: float
    seed_error_m: float
    ray: RayTrace
    max_path_alt_km: float

    def to_dict(self) -> dict:
        return {
            "label": self.label,
            "ox_mode": self.ox_mode,
            "frequency_mhz": self.frequency_mhz,
            "seed_error_m": self.seed_error_m,
            **{key: value for key, value in self.ray.to_dict().items() if key != "frequency_mhz"},
            "max_path_alt_km": self.max_path_alt_km,
        }


@dataclass(frozen=True)
class ModeSearchResult:
    best: ModeSolution
    homed: tuple[ModeSolution, ...]

    def to_dict(self) -> dict:
        return {
            "best": self.best.to_dict(),
            "homed": [solution.to_dict() for solution in self.homed],
        }


@dataclass(frozen=True)
class FrequencySweepResult:
    frequency_mhz: float
    ordinary: ModeSolution
    extraordinary: ModeSolution
    ordinary_homed: tuple[ModeSolution, ...] = ()
    extraordinary_homed: tuple[ModeSolution, ...] = ()

    def to_dict(self) -> dict:
        return {
            "frequency_mhz": self.frequency_mhz,
            "o_mode": self.ordinary.to_dict(),
            "x_mode": self.extraordinary.to_dict(),
            "o_mode_homed": [solution.to_dict() for solution in self.ordinary_homed],
            "x_mode_homed": [solution.to_dict() for solution in self.extraordinary_homed],
        }


@dataclass(frozen=True)
class SpaceToSpaceDemoResult:
    scenario: SpaceToSpaceScenario
    rx: GeoPoint
    line_of_sight_bearing_deg: float
    line_of_sight_elevation_deg: float
    line_of_sight_slant_range_km: float
    sweep: tuple[FrequencySweepResult, ...]

    def to_dict(self) -> dict:
        return {
            "when": self.scenario.when.isoformat(),
            "tx": {
                "lat_deg": self.scenario.tx.lat_deg,
                "lon_deg": self.scenario.tx.lon_deg,
                "alt_km": self.scenario.tx.alt_km,
            },
            "rx": {
                "lat_deg": self.rx.lat_deg,
                "lon_deg": self.rx.lon_deg,
                "alt_km": self.rx.alt_km,
            },
            "link": {
                "distance_km": self.scenario.distance_km,
                "bearing_deg": self.scenario.bearing_deg,
                "frequency_start_mhz": self.scenario.frequency_start_mhz,
                "frequency_stop_mhz": self.scenario.frequency_stop_mhz,
                "frequency_step_khz": self.scenario.frequency_step_khz,
                "count": len(self.sweep),
                "nhops": self.scenario.nhops,
            },
            "line_of_sight": {
                "bearing_deg": self.line_of_sight_bearing_deg,
                "elevation_deg": self.line_of_sight_elevation_deg,
                "slant_range_km": self.line_of_sight_slant_range_km,
            },
            "sweep": [entry.to_dict() for entry in self.sweep],
        }


def default_scenario() -> SpaceToSpaceScenario:
    return SpaceToSpaceScenario(
        when=dt.datetime(2020, 1, 15, 0, 0, 0),
        tx=GeoPoint(34.873, -106.614, 500.0),
    )


def _parse_time(text: str) -> dt.datetime:
    when = dt.datetime.fromisoformat(text)
    if when.tzinfo is not None:
        when = when.astimezone(dt.timezone.utc).replace(tzinfo=None)
    return when


def _line_of_sight_angles(tx: GeoPoint, rx: GeoPoint) -> tuple[float, float, float]:
    tx_xyz = np.asarray(llh_to_ecef(tx.lat_deg, tx.lon_deg, tx.alt_km * 1000.0), dtype=float)
    rx_xyz = np.asarray(llh_to_ecef(rx.lat_deg, rx.lon_deg, rx.alt_km * 1000.0), dtype=float)
    if tx_xyz.ndim > 1:
        tx_xyz = tx_xyz[0]
    if rx_xyz.ndim > 1:
        rx_xyz = rx_xyz[0]
    los = rx_xyz - tx_xyz
    slant_range_km = float(np.linalg.norm(los) / 1000.0)

    lat_rad = math.radians(tx.lat_deg)
    lon_rad = math.radians(tx.lon_deg)
    east = np.array([-math.sin(lon_rad), math.cos(lon_rad), 0.0], dtype=float)
    north = np.array(
        [
            -math.sin(lat_rad) * math.cos(lon_rad),
            -math.sin(lat_rad) * math.sin(lon_rad),
            math.cos(lat_rad),
        ],
        dtype=float,
    )
    up = np.array(
        [
            math.cos(lat_rad) * math.cos(lon_rad),
            math.cos(lat_rad) * math.sin(lon_rad),
            math.sin(lat_rad),
        ],
        dtype=float,
    )
    east_component = float(np.dot(los, east))
    north_component = float(np.dot(los, north))
    up_component = float(np.dot(los, up))
    bearing_deg = math.degrees(math.atan2(east_component, north_component))
    elevation_deg = math.degrees(math.atan2(up_component, math.hypot(east_component, north_component)))
    return bearing_deg, elevation_deg, slant_range_km


def _raytrace_grid_kwargs() -> dict[str, float | int | str]:
    return {
        "alt_min_km": 60.0,
        "alt_max_km": 900.0,
        "alt_step_km": 5.0,
        "lat_step_deg": 0.5,
        "lon_step_deg": 0.5,
        "lat_margin_deg": 3.0,
        "lon_margin_deg": 3.0,
        "d_region_model": "fpt2018",
    }


def _grid_cache_key(scenario: SpaceToSpaceScenario, rx: GeoPoint) -> str:
    payload = {
        "when": scenario.when.isoformat(),
        "tx": [scenario.tx.lat_deg, scenario.tx.lon_deg, scenario.tx.alt_km],
        "rx": [rx.lat_deg, rx.lon_deg, rx.alt_km],
        "f107": scenario.f107,
        "ap_daily": scenario.ap_daily,
        "grid": _raytrace_grid_kwargs(),
    }
    return hashlib.sha1(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]


def _default_grid_cache_path(scenario: SpaceToSpaceScenario, rx: GeoPoint) -> Path:
    return Path(".cache") / "space_to_space_grids" / f"{_grid_cache_key(scenario, rx)}.npz"


def _load_or_build_demo_grid(
    scenario: SpaceToSpaceScenario,
    rx: GeoPoint,
    *,
    grid_cache_path: Path | None,
    rebuild_grid: bool,
) -> IonosphereGrid:
    cache_path = _default_grid_cache_path(scenario, rx) if grid_cache_path is None else grid_cache_path.expanduser()
    if cache_path.exists() and not rebuild_grid:
        return load_ionosphere_grid(cache_path)
    grid = build_pyiri_grid(
        scenario.when,
        scenario.tx,
        rx,
        f107=scenario.f107,
        ap_daily=scenario.ap_daily,
        **_raytrace_grid_kwargs(),
    )
    save_ionosphere_grid(cache_path, grid)
    return grid


def _surface_distance_km(start: GeoPoint, end: GeoPoint) -> float:
    lat1 = math.radians(start.lat_deg)
    lat2 = math.radians(end.lat_deg)
    dlat = lat2 - lat1
    dlon = math.radians(((end.lon_deg - start.lon_deg + 180.0) % 360.0) - 180.0)
    a = math.sin(dlat / 2.0) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2
    return 6371.0088 * 2.0 * math.asin(min(1.0, math.sqrt(max(a, 0.0))))


def _max_path_alt_km(ray: RayTrace) -> float:
    heights = np.asarray(ray.path.get("height", []), dtype=float)
    heights = heights[np.isfinite(heights) & (heights < 1e40)]
    return float(np.max(heights)) if heights.size else float("nan")


def _trace_explicit_state_vector(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid,
    elevation_deg: float,
    bearing_deg: float,
    frequency_mhz: float,
    ox_mode: int,
    nhops: int,
) -> RayTrace:
    prepared = tracer.prepare_ray_state_vector_batch(
        tx=tx,
        grid=grid,
        elevations_deg=[elevation_deg],
        bearings_deg=[bearing_deg],
        freqs_mhz=[frequency_mhz],
        ox_mode=ox_mode,
    )
    if prepared is None:
        return RayTrace(summary={}, path={"initial_elev": math.nan, "initial_bearing": math.nan, "frequency": frequency_mhz}, state={})
    tx_local, state_vector, valid_mask = prepared
    if valid_mask.size != 1 or not bool(valid_mask[0]):
        return RayTrace(summary={}, path={"initial_elev": math.nan, "initial_bearing": math.nan, "frequency": frequency_mhz}, state={})
    ray = tracer.trace_state_vector_batch(
        tx=tx_local,
        elevations_deg=[elevation_deg],
        bearings_deg=[bearing_deg],
        freqs_mhz=[frequency_mhz],
        grid=grid,
        state_vector=state_vector,
        ox_mode=ox_mode,
        nhops=nhops,
    )[0]
    distance = ray_point_distance(ray.path, rx, expect_reflection=False)
    ray.error_m = distance.distance_m
    if math.isfinite(distance.distance_m):
        ray.group_range_to_rx_km = distance.group_path_km
        ray.geometric_dist_to_rx_km = distance.geometric_path_km
        ray.total_absorption_db = distance.total_absorption_db
    return ray


def _trace_explicit_state_vector_batch(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid,
    elevations_deg: np.ndarray,
    bearings_deg: np.ndarray,
    frequency_mhz: float,
    ox_mode: int,
    nhops: int,
) -> list[RayTrace]:
    elevs = np.asarray(elevations_deg, dtype=float)
    bears = np.asarray(bearings_deg, dtype=float)
    freqs = np.full(elevs.shape, float(frequency_mhz), dtype=float)
    prepared = tracer.prepare_ray_state_vector_batch(
        tx=tx,
        grid=grid,
        elevations_deg=elevs,
        bearings_deg=bears,
        freqs_mhz=freqs,
        ox_mode=ox_mode,
    )
    if prepared is None:
        return []
    tx_local, state_vector, valid_mask = prepared
    rays = tracer.trace_state_vector_batch(
        tx=tx_local,
        elevations_deg=elevs[valid_mask],
        bearings_deg=bears[valid_mask],
        freqs_mhz=freqs[valid_mask],
        grid=grid,
        state_vector=state_vector,
        ox_mode=ox_mode,
        nhops=nhops,
    )
    output: list[RayTrace] = []
    ray_iter = iter(rays)
    for is_valid, elev_deg, bear_deg in zip(valid_mask, elevs, bears):
        if not bool(is_valid):
            output.append(
                RayTrace(
                    summary={},
                    path={"initial_elev": float(elev_deg), "initial_bearing": float(bear_deg), "frequency": frequency_mhz},
                    state={},
                )
            )
            continue
        ray = next(ray_iter)
        distance = ray_point_distance(ray.path, rx, expect_reflection=False)
        ray.error_m = distance.distance_m
        if math.isfinite(distance.distance_m):
            ray.group_range_to_rx_km = distance.group_path_km
            ray.geometric_dist_to_rx_km = distance.geometric_path_km
            ray.total_absorption_db = distance.total_absorption_db
        output.append(ray)
    return output


def _mark_homed(ray: RayTrace, homing_tolerance_m: float) -> RayTrace:
    if ray.error_m < homing_tolerance_m:
        ray.home = True
        ray.perigee_km = float(np.nanmin(np.asarray(ray.path.get("height", []), dtype=float)))
    return ray


def _coerce_lon_to_target(lon_deg: float, target_lon_deg: float) -> float:
    candidates = np.array([lon_deg - 360.0, lon_deg, lon_deg + 360.0], dtype=float)
    return float(candidates[np.argmin(np.abs(candidates - target_lon_deg))])


def _ray_shell_pierce_point(ray: RayTrace, shell_alt_km: float) -> tuple[float, float] | None:
    heights = np.asarray(ray.path.get("height", []), dtype=float)
    lats = np.asarray(ray.path.get("lat", []), dtype=float)
    lons = np.asarray(ray.path.get("lon", []), dtype=float)
    valid = np.isfinite(heights) & np.isfinite(lats) & np.isfinite(lons) & (heights < 1e40)
    heights = heights[valid]
    lats = lats[valid]
    lons = lons[valid]
    if heights.size < 2:
        return None
    crossing: tuple[float, float] | None = None
    for idx in range(heights.size - 1):
        h0 = float(heights[idx])
        h1 = float(heights[idx + 1])
        if not (h1 > h0):
            continue
        if shell_alt_km < h0 or shell_alt_km > h1:
            continue
        t = 0.0 if h1 == h0 else float((shell_alt_km - h0) / (h1 - h0))
        lat = float(lats[idx] + t * (lats[idx + 1] - lats[idx]))
        lon1 = _coerce_lon_to_target(float(lons[idx + 1]), float(lons[idx]))
        lon = float(lons[idx] + t * (lon1 - float(lons[idx])))
        crossing = (lat, lon)
    return crossing


def _fan_interpolated_update(
    rays: list[RayTrace],
    rx: GeoPoint,
) -> tuple[float, float] | None:
    points: list[tuple[float, float]] = []
    launch_bears: list[float] = []
    launch_elevs: list[float] = []
    for ray in rays:
        pierce = _ray_shell_pierce_point(ray, rx.alt_km)
        if pierce is None:
            continue
        points.append((pierce[0], _coerce_lon_to_target(pierce[1], rx.lon_deg)))
        launch_bears.append(float(ray.launch_bearing_deg))
        launch_elevs.append(float(ray.launch_elevation_deg))
    if len(points) < 3:
        return None
    point_array = np.asarray(points, dtype=float)
    rx_lon = _coerce_lon_to_target(rx.lon_deg, float(np.mean(point_array[:, 1])))
    target = np.array([[rx.lat_deg, rx_lon]], dtype=float)
    bearing = griddata(point_array, np.asarray(launch_bears, dtype=float), target, method="linear")
    elevation = griddata(point_array, np.asarray(launch_elevs, dtype=float), target, method="linear")
    if bearing is None or elevation is None or not np.isfinite(bearing[0]) or not np.isfinite(elevation[0]):
        bearing = griddata(point_array, np.asarray(launch_bears, dtype=float), target, method="nearest")
        elevation = griddata(point_array, np.asarray(launch_elevs, dtype=float), target, method="nearest")
    if bearing is None or elevation is None or not np.isfinite(bearing[0]) or not np.isfinite(elevation[0]):
        return None
    return float(elevation[0]), float(bearing[0])


def _pattern_refine_branch(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid,
    start_ray: RayTrace,
    frequency_mhz: float,
    ox_mode: int,
    nhops: int,
    homing_tolerance_m: float,
    ray_cache: dict[tuple[int, int], RayTrace],
    max_iter: int = 10,
) -> RayTrace:
    step_elev_deg = 0.5
    step_bear_deg = 0.5
    center_elev = float(start_ray.launch_elevation_deg)
    center_bear = float(start_ray.launch_bearing_deg)
    best = start_ray

    def cache_key(elev_deg: float, bear_deg: float) -> tuple[int, int]:
        return (int(round(elev_deg * 1e5)), int(round(bear_deg * 1e5)))

    for _ in range(max_iter):
        elev_offsets = np.array([-step_elev_deg, 0.0, step_elev_deg], dtype=float)
        bear_offsets = np.array([-step_bear_deg, 0.0, step_bear_deg], dtype=float)
        elev_grid, bear_grid = np.meshgrid(center_elev + elev_offsets, center_bear + bear_offsets, indexing="xy")
        candidate_elevs = elev_grid.ravel()
        candidate_bears = bear_grid.ravel()

        uncached_mask = np.array([cache_key(elev, bear) not in ray_cache for elev, bear in zip(candidate_elevs, candidate_bears)], dtype=bool)
        if np.any(uncached_mask):
            new_rays = _trace_explicit_state_vector_batch(
                tracer,
                tx=tx,
                rx=rx,
                grid=grid,
                elevations_deg=candidate_elevs[uncached_mask],
                bearings_deg=candidate_bears[uncached_mask],
                frequency_mhz=frequency_mhz,
                ox_mode=ox_mode,
                nhops=nhops,
            )
            for elev_deg, bear_deg, ray in zip(candidate_elevs[uncached_mask], candidate_bears[uncached_mask], new_rays):
                ray_cache[cache_key(float(elev_deg), float(bear_deg))] = ray

        rays = [ray_cache[cache_key(float(elev_deg), float(bear_deg))] for elev_deg, bear_deg in zip(candidate_elevs, candidate_bears)]
        finite_rays = [ray for ray in rays if math.isfinite(ray.error_m)]
        if not finite_rays:
            break
        iter_best = min(finite_rays, key=lambda ray: ray.error_m)
        if iter_best.error_m < best.error_m:
            best = iter_best
        if best.error_m < homing_tolerance_m:
            return _mark_homed(best, homing_tolerance_m)
        updated = _fan_interpolated_update(finite_rays, rx)
        if updated is None:
            center_elev = float(iter_best.launch_elevation_deg)
            center_bear = float(iter_best.launch_bearing_deg)
            step_elev_deg *= 0.7
            step_bear_deg *= 0.7
        else:
            center_elev, center_bear = updated
            best_key = cache_key(center_elev, center_bear)
            cached = ray_cache.get(best_key)
            if cached is not None and cached.error_m < best.error_m:
                best = cached
        if step_elev_deg < 0.02 and step_bear_deg < 0.02:
            break
    return _mark_homed(best, homing_tolerance_m)


def _solve_mode(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid,
    frequency_mhz: float,
    ox_mode: int,
    nhops: int,
    homing_tolerance_m: float,
    start_angles_deg: tuple[float, float],
) -> ModeSearchResult:
    start_candidates = _build_start_candidates(
        tracer,
        tx=tx,
        rx=rx,
        grid=grid,
        frequency_mhz=frequency_mhz,
        ox_mode=ox_mode,
        nhops=nhops,
        start_angles_deg=start_angles_deg,
    )
    if not start_candidates:
        empty = ModeSolution(
            label="O-mode" if ox_mode == 1 else "X-mode",
            ox_mode=ox_mode,
            frequency_mhz=frequency_mhz,
            seed_error_m=math.inf,
            ray=RayTrace(summary={}, path={"initial_elev": math.nan, "initial_bearing": math.nan, "frequency": frequency_mhz}, state={}),
            max_path_alt_km=float("nan"),
        )
        return ModeSearchResult(best=empty, homed=())

    ray_cache: dict[tuple[int, int], RayTrace] = {}

    def cached_trace(elevation_deg: float, bearing_deg: float) -> RayTrace:
        key = (int(round(elevation_deg * 1e5)), int(round(bearing_deg * 1e5)))
        cached = ray_cache.get(key)
        if cached is not None:
            return cached
        ray = _trace_explicit_state_vector(
            tracer,
            tx=tx,
            rx=rx,
            grid=grid,
            elevation_deg=elevation_deg,
            bearing_deg=bearing_deg,
            frequency_mhz=frequency_mhz,
            ox_mode=ox_mode,
            nhops=nhops,
        )
        ray_cache[key] = ray
        return ray

    attempts: list[ModeSolution] = []
    label = "O-mode" if ox_mode == 1 else "X-mode"
    for start_ray in start_candidates:
        if start_ray.error_m < homing_tolerance_m:
            ray = _mark_homed(start_ray, homing_tolerance_m)
        else:
            ray = _pattern_refine_branch(
                tracer,
                tx=tx,
                rx=rx,
                grid=grid,
                start_ray=start_ray,
                frequency_mhz=frequency_mhz,
                ox_mode=ox_mode,
                nhops=nhops,
                homing_tolerance_m=homing_tolerance_m,
                ray_cache=ray_cache,
            )
        if ray.error_m < homing_tolerance_m:
            ray = _mark_homed(ray, homing_tolerance_m)
        attempts.append(
            ModeSolution(
                label=label,
                ox_mode=ox_mode,
                frequency_mhz=frequency_mhz,
                seed_error_m=start_ray.error_m,
                ray=ray,
                max_path_alt_km=_max_path_alt_km(ray),
            )
        )
    best = min(attempts, key=lambda item: item.ray.error_m)
    homed = _dedupe_homed_solutions([attempt for attempt in attempts if attempt.ray.home])
    return ModeSearchResult(
        best=best,
        homed=tuple(homed),
    )


def _topside_search_arrays(tx: GeoPoint, rx: GeoPoint) -> tuple[np.ndarray, np.ndarray]:
    gc_bearing_deg, gc_elevation_deg, slant_range_km = _line_of_sight_angles(tx, rx)
    tx_ground = GeoPoint(tx.lat_deg, tx.lon_deg, 0.0)
    rx_ground = GeoPoint(rx.lat_deg, rx.lon_deg, 0.0)
    mid_lats, mid_lons = great_circle_waypoints(tx_ground, rx_ground, count=3)
    midpoint_ground = GeoPoint(float(mid_lats[1]), float(mid_lons[1]), 0.0)
    tangent_elevation_deg = _line_of_sight_angles(tx, midpoint_ground)[1]
    if not math.isfinite(tangent_elevation_deg):
        tangent_elevation_deg = gc_elevation_deg - 40.0
    elevation_min = max(-89.0, min(tangent_elevation_deg, gc_elevation_deg))
    elevation_max = min(89.0, max(tangent_elevation_deg, gc_elevation_deg))
    elevations_deg = np.linspace(elevation_min, elevation_max, 41, dtype=float)
    az_half_width_deg = max(3.0, float(math.ceil(30.0 * math.sqrt(100.0 / max(slant_range_km, 100.0)))))
    bearings_deg = np.linspace(gc_bearing_deg - az_half_width_deg, gc_bearing_deg + az_half_width_deg, 17, dtype=float)
    return elevations_deg, bearings_deg


def _coarse_search_candidates(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid,
    frequency_mhz: float,
    ox_mode: int,
    nhops: int,
) -> list[RayTrace]:
    elevations_deg, bearings_deg = _topside_search_arrays(tx, rx)
    elevation_grid, bearing_grid = np.meshgrid(elevations_deg, bearings_deg, indexing="xy")
    elevs = elevation_grid.ravel()
    bears = bearing_grid.ravel()
    freqs = np.full(elevs.shape, float(frequency_mhz), dtype=float)
    prepared = tracer.prepare_ray_state_vector_batch(
        tx=tx,
        grid=grid,
        elevations_deg=elevs,
        bearings_deg=bears,
        freqs_mhz=freqs,
        ox_mode=ox_mode,
    )
    if prepared is None:
        return []
    tx_local, state_vector, valid_mask = prepared
    if not np.any(valid_mask):
        return []
    rays = tracer.trace_state_vector_batch(
        tx=tx_local,
        elevations_deg=elevs[valid_mask],
        bearings_deg=bears[valid_mask],
        freqs_mhz=freqs[valid_mask],
        grid=grid,
        state_vector=state_vector,
        ox_mode=ox_mode,
        nhops=nhops,
    )
    valid_rays: list[RayTrace] = []
    for ray in rays:
        distance = ray_point_distance(ray.path, rx, expect_reflection=False)
        if not math.isfinite(distance.distance_m):
            continue
        ray.error_m = distance.distance_m
        ray.group_range_to_rx_km = distance.group_path_km
        ray.geometric_dist_to_rx_km = distance.geometric_path_km
        ray.total_absorption_db = distance.total_absorption_db
        valid_rays.append(ray)
    valid_rays.sort(key=lambda ray: ray.error_m)
    return valid_rays


def _dedupe_start_candidates(rays: list[RayTrace], max_count: int = 6) -> list[RayTrace]:
    chosen: list[RayTrace] = []
    seen: set[tuple[int, int]] = set()
    for ray in sorted(rays, key=lambda item: item.error_m):
        elev = ray.launch_elevation_deg
        bearing = ray.launch_bearing_deg
        if not (math.isfinite(elev) and math.isfinite(bearing)):
            continue
        key = (int(round(elev * 10.0)), int(round(bearing * 10.0)))
        if key in seen:
            continue
        seen.add(key)
        chosen.append(ray)
        if len(chosen) >= max_count:
            break
    return chosen


def _local_minimum_mask(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    if values.ndim != 2 or not np.any(finite):
        return np.zeros_like(values, dtype=bool)
    padded = np.pad(values, 1, mode="constant", constant_values=np.inf)
    center = padded[1:-1, 1:-1]
    is_local = finite.copy()
    for di in range(3):
        for dj in range(3):
            if di == 1 and dj == 1:
                continue
            neighbor = padded[di:di + values.shape[0], dj:dj + values.shape[1]]
            is_local &= center <= neighbor
    return is_local


def _coarse_local_minima_candidates(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid,
    frequency_mhz: float,
    ox_mode: int,
    nhops: int,
    max_count: int = 8,
) -> list[RayTrace]:
    elevations_deg, bearings_deg = _topside_search_arrays(tx, rx)
    elevation_grid, bearing_grid = np.meshgrid(elevations_deg, bearings_deg, indexing="xy")
    elevs = elevation_grid.ravel()
    bears = bearing_grid.ravel()
    freqs = np.full(elevs.shape, float(frequency_mhz), dtype=float)
    prepared = tracer.prepare_ray_state_vector_batch(
        tx=tx,
        grid=grid,
        elevations_deg=elevs,
        bearings_deg=bears,
        freqs_mhz=freqs,
        ox_mode=ox_mode,
    )
    if prepared is None:
        return []
    tx_local, state_vector, valid_mask = prepared
    if not np.any(valid_mask):
        return []
    traced = tracer.trace_state_vector_batch(
        tx=tx_local,
        elevations_deg=elevs[valid_mask],
        bearings_deg=bears[valid_mask],
        freqs_mhz=freqs[valid_mask],
        grid=grid,
        state_vector=state_vector,
        ox_mode=ox_mode,
        nhops=nhops,
    )
    error_grid = np.full(elevation_grid.shape, np.inf, dtype=float)
    valid_rays: list[RayTrace] = []
    ray_lookup: dict[tuple[int, int], RayTrace] = {}
    for flat_index, ray in zip(np.flatnonzero(valid_mask), traced):
        distance = ray_point_distance(ray.path, rx, expect_reflection=False)
        if not math.isfinite(distance.distance_m):
            continue
        ray.error_m = distance.distance_m
        ray.group_range_to_rx_km = distance.group_path_km
        ray.geometric_dist_to_rx_km = distance.geometric_path_km
        ray.total_absorption_db = distance.total_absorption_db
        ij = np.unravel_index(int(flat_index), elevation_grid.shape)
        error_grid[ij] = distance.distance_m
        ray_lookup[ij] = ray
        valid_rays.append(ray)
    if not valid_rays:
        return []
    local_mask = _local_minimum_mask(error_grid)
    minima = [ray_lookup[index] for index in ray_lookup if local_mask[index]]
    if not minima:
        minima = _select_family_candidates(valid_rays, max_count=max_count)
    return _dedupe_start_candidates(sorted(minima, key=lambda item: item.error_m), max_count=max_count)


def _select_family_candidates(rays: list[RayTrace], max_count: int = 8) -> list[RayTrace]:
    if not rays:
        return []
    family_best: dict[int, RayTrace] = {}
    for ray in rays:
        elev = ray.launch_elevation_deg
        bearing = ray.launch_bearing_deg
        if not (math.isfinite(elev) and math.isfinite(bearing)):
            continue
        family_key = int(math.floor((elev + 90.0) / 4.0))
        incumbent = family_best.get(family_key)
        if incumbent is None or ray.error_m < incumbent.error_m:
            family_best[family_key] = ray
    selected = sorted(family_best.values(), key=lambda item: item.error_m)
    if selected:
        global_best = min(rays, key=lambda item: item.error_m)
        if all(global_best is not item for item in selected):
            selected.insert(0, global_best)
    return _dedupe_start_candidates(selected, max_count=max_count)


def _dedupe_homed_solutions(solutions: list[ModeSolution], max_count: int = 12) -> list[ModeSolution]:
    chosen: list[ModeSolution] = []
    seen: set[tuple[int, int]] = set()
    for solution in sorted(
        solutions,
        key=lambda item: (
            item.ray.group_range_to_rx_km if item.ray.group_range_to_rx_km is not None else math.inf,
            item.ray.error_m,
        ),
    ):
        elev = solution.ray.launch_elevation_deg
        group_range = solution.ray.group_range_to_rx_km
        if not (math.isfinite(elev) and group_range is not None and math.isfinite(group_range)):
            continue
        key = (
            int(round(elev * 5.0)),
            int(round(group_range)),
        )
        if key in seen:
            continue
        seen.add(key)
        chosen.append(solution)
        if len(chosen) >= max_count:
            break
    return chosen


def _build_start_candidates(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid,
    frequency_mhz: float,
    ox_mode: int,
    nhops: int,
    start_angles_deg: tuple[float, float],
) -> list[RayTrace]:
    candidates: list[RayTrace] = []
    seeded = _trace_explicit_state_vector(
        tracer,
        tx=tx,
        rx=rx,
        grid=grid,
        elevation_deg=float(start_angles_deg[0]),
        bearing_deg=float(start_angles_deg[1]),
        frequency_mhz=frequency_mhz,
        ox_mode=ox_mode,
        nhops=nhops,
    )
    if math.isfinite(seeded.error_m):
        candidates.append(seeded)
    coarse_candidates = _coarse_local_minima_candidates(
        tracer,
        tx=tx,
        rx=rx,
        grid=grid,
        frequency_mhz=frequency_mhz,
        ox_mode=ox_mode,
        nhops=nhops,
        max_count=8,
    )
    candidates.extend(coarse_candidates)
    return _dedupe_start_candidates(candidates)


def _next_seed(result: ModeSearchResult, fallback: tuple[float, float]) -> tuple[float, float]:
    solution = result.best
    elev = solution.ray.launch_elevation_deg
    bearing = solution.ray.launch_bearing_deg
    if not (math.isfinite(elev) and math.isfinite(bearing)):
        return fallback
    if solution.ray.home or solution.ray.error_m <= solution.seed_error_m:
        return float(elev), float(bearing)
    return fallback


def _solve_frequency_sweep(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid,
    frequencies_mhz: np.ndarray,
    nhops: int,
    homing_tolerance_m: float,
    start_angles_deg: tuple[float, float],
) -> tuple[FrequencySweepResult, ...]:
    seeds = {
        1: tuple(float(value) for value in start_angles_deg),
        -1: tuple(float(value) for value in start_angles_deg),
    }
    solved_descending: list[FrequencySweepResult] = []
    for frequency_mhz in sorted((float(value) for value in frequencies_mhz), reverse=True):
        ordinary = _solve_mode(
            tracer,
            tx=tx,
            rx=rx,
            grid=grid,
            frequency_mhz=frequency_mhz,
            ox_mode=1,
            nhops=nhops,
            homing_tolerance_m=homing_tolerance_m,
            start_angles_deg=seeds[1],
        )
        seeds[1] = _next_seed(ordinary, seeds[1])
        extraordinary = _solve_mode(
            tracer,
            tx=tx,
            rx=rx,
            grid=grid,
            frequency_mhz=frequency_mhz,
            ox_mode=-1,
            nhops=nhops,
            homing_tolerance_m=homing_tolerance_m,
            start_angles_deg=seeds[-1],
        )
        seeds[-1] = _next_seed(extraordinary, seeds[-1])
        solved_descending.append(
            FrequencySweepResult(
                frequency_mhz=frequency_mhz,
                ordinary=ordinary.best,
                extraordinary=extraordinary.best,
                ordinary_homed=ordinary.homed,
                extraordinary_homed=extraordinary.homed,
            )
        )
    return tuple(sorted(solved_descending, key=lambda entry: entry.frequency_mhz))


def _ray_profile(
    tx: GeoPoint,
    ray: RayTrace,
    *,
    longitude_reference_deg: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    lats = np.asarray(ray.path.get("lat", []), dtype=float)
    lons = np.asarray(ray.path.get("lon", []), dtype=float)
    alts = np.asarray(ray.path.get("height", []), dtype=float)
    valid = np.isfinite(lats) & np.isfinite(lons) & np.isfinite(alts) & (alts < 1e40)
    lats = lats[valid]
    if longitude_reference_deg is None:
        lons = np.asarray([wrap_longitude(float(value)) for value in lons[valid]], dtype=float)
    else:
        lons = np.asarray(
            [coerce_longitude_for_grid(float(value), longitude_reference_deg) for value in lons[valid]],
            dtype=float,
        )
    alts = alts[valid]
    along_km = np.asarray(
        [_surface_distance_km(tx, GeoPoint(float(lat), float(lon), 0.0)) for lat, lon in zip(lats, lons)],
        dtype=float,
    )
    return lats, lons, np.maximum.accumulate(along_km), alts


def _centers_to_edges(values: np.ndarray) -> np.ndarray:
    centers = np.asarray(values, dtype=float)
    if centers.size == 1:
        return np.array([centers[0] - 0.5, centers[0] + 0.5], dtype=float)
    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def _visible_color_limits(
    field: np.ndarray,
    *,
    x_centers: np.ndarray,
    y_centers: np.ndarray,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
) -> tuple[float, float]:
    xmin, xmax = sorted((float(x_limits[0]), float(x_limits[1])))
    ymin, ymax = sorted((float(y_limits[0]), float(y_limits[1])))
    xmask = (np.asarray(x_centers, dtype=float) >= xmin) & (np.asarray(x_centers, dtype=float) <= xmax)
    ymask = (np.asarray(y_centers, dtype=float) >= ymin) & (np.asarray(y_centers, dtype=float) <= ymax)
    if not np.any(xmask) or not np.any(ymask):
        finite = np.asarray(field, dtype=float)[np.isfinite(field)]
    else:
        finite = np.asarray(field, dtype=float)[np.ix_(ymask, xmask)]
        finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 0.0, 1.0
    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    if not vmax > vmin:
        vmax = vmin + 1e-6
    return vmin, vmax


def _homed_mode_solutions(result: SpaceToSpaceDemoResult, ox_mode: int) -> tuple[ModeSolution, ...]:
    if ox_mode == 1:
        return tuple(solution for entry in result.sweep for solution in entry.ordinary_homed)
    return tuple(solution for entry in result.sweep for solution in entry.extraordinary_homed)


def _plan_extent(
    result: SpaceToSpaceDemoResult,
    *,
    grid_longitudes_deg: np.ndarray,
    homed_modes: tuple[ModeSolution, ...],
) -> tuple[tuple[float, float], tuple[float, float]]:
    lat_values = [result.scenario.tx.lat_deg, result.rx.lat_deg]
    lon_values = [
        coerce_longitude_for_grid(result.scenario.tx.lon_deg, grid_longitudes_deg),
        coerce_longitude_for_grid(result.rx.lon_deg, grid_longitudes_deg),
    ]
    for solution in homed_modes:
        lats, lons, _, _ = _ray_profile(result.scenario.tx, solution.ray, longitude_reference_deg=grid_longitudes_deg)
        if lats.size:
            lat_values.extend(float(value) for value in lats)
        if lons.size:
            lon_values.extend(float(value) for value in lons)
    lat_min = min(lat_values)
    lat_max = max(lat_values)
    lon_min = min(lon_values)
    lon_max = max(lon_values)
    lat_margin = max(0.2, 0.08 * max(lat_max - lat_min, 0.5))
    lon_margin = max(0.2, 0.08 * max(lon_max - lon_min, 0.5))
    return (lon_min - lon_margin, lon_max + lon_margin), (lat_min - lat_margin, lat_max + lat_margin)


def _ionogram_series(mode_solutions: tuple[ModeSolution, ...]) -> tuple[np.ndarray, np.ndarray]:
    freqs: list[float] = []
    virtual_ranges: list[float] = []
    for solution in mode_solutions:
        group_range_km = solution.ray.group_range_to_rx_km
        if group_range_km is None or not math.isfinite(group_range_km):
            continue
        freqs.append(float(solution.frequency_mhz))
        virtual_ranges.append(float(group_range_km))
    return np.asarray(freqs, dtype=float), np.asarray(virtual_ranges, dtype=float)


def plot_space_to_space_overview(
    path: Path,
    result: SpaceToSpaceDemoResult,
    *,
    grid_cache_path: Path | None = None,
    rebuild_grid: bool = False,
) -> Path:
    import matplotlib.pyplot as plt  # type: ignore
    from matplotlib.colors import Normalize  # type: ignore
    from matplotlib.lines import Line2D  # type: ignore

    grid = _load_or_build_demo_grid(
        result.scenario,
        result.rx,
        grid_cache_path=grid_cache_path,
        rebuild_grid=rebuild_grid,
    )
    homed_o = _homed_mode_solutions(result, 1)
    homed_x = _homed_mode_solutions(result, -1)
    homed_all = homed_o + homed_x
    density_m3 = np.asarray(grid.iono_en_grid, dtype=float) * 1e6
    peak_density_m3 = np.max(density_m3, axis=2)
    fof2_mhz = 8.98e-6 * np.sqrt(np.maximum(peak_density_m3, 0.0))
    lon_edges = _centers_to_edges(grid.longitudes_deg)
    lat_edges = _centers_to_edges(grid.latitudes_deg)
    lon_mesh, lat_mesh = np.meshgrid(lon_edges, lat_edges, indexing="xy")
    tx_lon_plot = coerce_longitude_for_grid(result.scenario.tx.lon_deg, grid.longitudes_deg)
    rx_lon_plot = coerce_longitude_for_grid(result.rx.lon_deg, grid.longitudes_deg)

    tx_ground = GeoPoint(result.scenario.tx.lat_deg, result.scenario.tx.lon_deg, 0.0)
    rx_ground = GeoPoint(result.rx.lat_deg, result.rx.lon_deg, 0.0)
    track_lats_deg, track_lons_deg = great_circle_waypoints(tx_ground, rx_ground, count=256)
    along_track_km = np.linspace(0.0, result.scenario.distance_km, track_lats_deg.size, dtype=float)
    sample_lons_deg = np.asarray(
        [coerce_longitude_for_grid(float(lon), grid.longitudes_deg) for lon in track_lons_deg],
        dtype=float,
    )
    lat_mesh_track, alt_mesh_track = np.meshgrid(np.asarray(track_lats_deg, dtype=float), grid.altitudes_km, indexing="xy")
    lon_mesh_track, _ = np.meshgrid(sample_lons_deg, grid.altitudes_km, indexing="xy")
    swath_points = np.column_stack((lat_mesh_track.ravel(), lon_mesh_track.ravel(), alt_mesh_track.ravel()))
    density_interp = RegularGridInterpolator(
        (grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km),
        density_m3,
        bounds_error=False,
        fill_value=np.nan,
    )
    swath_density_m3 = density_interp(swath_points).reshape(grid.altitudes_km.size, along_track_km.size)
    swath_plasma_frequency_mhz = 8.98e-6 * np.sqrt(np.maximum(np.where(np.isfinite(swath_density_m3), swath_density_m3, 0.0), 0.0))
    along_edges = _centers_to_edges(along_track_km)
    alt_edges = _centers_to_edges(grid.altitudes_km)
    along_mesh, alt_mesh = np.meshgrid(along_edges, alt_edges, indexing="xy")

    xlim0, ylim0 = _plan_extent(result, grid_longitudes_deg=grid.longitudes_deg, homed_modes=homed_all)
    homed_max_altitudes = [solution.max_path_alt_km for solution in homed_all if math.isfinite(solution.max_path_alt_km)]
    max_alt_km = max(
        homed_max_altitudes + [result.scenario.tx.alt_km, result.rx.alt_km, float(grid.altitudes_km[0])],
    )
    xlim1 = (0.0, result.scenario.distance_km)
    ylim1 = (float(grid.altitudes_km[0]), max_alt_km + 20.0)
    fof2_vmin, fof2_vmax = _visible_color_limits(
        fof2_mhz,
        x_centers=np.asarray(grid.longitudes_deg, dtype=float),
        y_centers=np.asarray(grid.latitudes_deg, dtype=float),
        x_limits=xlim0,
        y_limits=ylim0,
    )
    pf_vmin, pf_vmax = _visible_color_limits(
        swath_plasma_frequency_mhz,
        x_centers=np.asarray(along_track_km, dtype=float),
        y_centers=np.asarray(grid.altitudes_km, dtype=float),
        x_limits=xlim1,
        y_limits=ylim1,
    )
    freq_norm = Normalize(vmin=float(result.scenario.frequency_start_mhz), vmax=float(result.scenario.frequency_stop_mhz))
    o_cmap = plt.get_cmap("Blues")
    x_cmap = plt.get_cmap("Oranges")
    iono_o_freqs, iono_o_ranges = _ionogram_series(homed_o)
    iono_x_freqs, iono_x_ranges = _ionogram_series(homed_x)
    ionogram_range_values = [result.line_of_sight_slant_range_km]
    if iono_o_ranges.size:
        ionogram_range_values.extend(float(value) for value in iono_o_ranges)
    if iono_x_ranges.size:
        ionogram_range_values.extend(float(value) for value in iono_x_ranges)
    iono_ymin = min(ionogram_range_values)
    iono_ymax = max(ionogram_range_values)
    iono_margin = max(10.0, 0.08 * max(iono_ymax - iono_ymin, 10.0))

    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, figsize=(18.2, 5.8))

    surface0 = ax0.pcolormesh(
        lon_mesh,
        lat_mesh,
        fof2_mhz,
        cmap="viridis",
        vmin=fof2_vmin,
        vmax=fof2_vmax,
        shading="flat",
        antialiased=False,
        rasterized=True,
    )
    for solution in homed_o:
        lats, lons, _, _ = _ray_profile(result.scenario.tx, solution.ray, longitude_reference_deg=grid.longitudes_deg)
        ax0.plot(lons, lats, color=o_cmap(0.35 + 0.6 * freq_norm(solution.frequency_mhz)), linewidth=1.8, alpha=0.95)
    for solution in homed_x:
        lats, lons, _, _ = _ray_profile(result.scenario.tx, solution.ray, longitude_reference_deg=grid.longitudes_deg)
        ax0.plot(lons, lats, color=x_cmap(0.35 + 0.6 * freq_norm(solution.frequency_mhz)), linewidth=1.7, alpha=0.92)
    ax0.scatter([tx_lon_plot], [result.scenario.tx.lat_deg], color="gold", edgecolors="black", s=170, marker="*", zorder=5, label="Tx")
    ax0.scatter([rx_lon_plot], [result.rx.lat_deg], color="tab:red", edgecolors="white", s=68, zorder=5, label="Rx")
    ax0.set_xlim(*xlim0)
    ax0.set_ylim(*ylim0)
    ax0.set_xlabel("Longitude (deg)")
    ax0.set_ylabel("Latitude (deg)")
    ax0.set_title("Plan View with foF2 / Peak Plasma Frequency")
    ax0.set_axisbelow(True)
    ax0.minorticks_on()
    ax0.grid(True, which="major", alpha=0.28, linewidth=0.7)
    ax0.grid(True, which="minor", alpha=0.12, linewidth=0.5, linestyle=":")
    ax0.legend(
        handles=[
            Line2D([0.0], [0.0], color=o_cmap(0.8), linewidth=2.0, label="O-mode homed"),
            Line2D([0.0], [0.0], color=x_cmap(0.8), linewidth=1.9, label="X-mode homed"),
            Line2D([0.0], [0.0], color="gold", marker="*", markersize=12, linewidth=0.0, markeredgecolor="black", label="Tx"),
            Line2D([0.0], [0.0], color="tab:red", marker="o", markersize=7, linewidth=0.0, markeredgecolor="white", label="Rx"),
        ],
        loc="best",
        framealpha=0.95,
    )
    cbar0 = fig.colorbar(surface0, ax=ax0, pad=0.02, label="foF2 / peak plasma frequency (MHz)")
    cbar0.ax.tick_params(labelsize=10)

    surface1 = ax1.pcolormesh(
        along_mesh,
        alt_mesh,
        swath_plasma_frequency_mhz,
        cmap="viridis",
        vmin=pf_vmin,
        vmax=pf_vmax,
        shading="flat",
        antialiased=False,
        rasterized=True,
    )
    for solution in homed_o:
        _, _, along_km, alts_km = _ray_profile(result.scenario.tx, solution.ray, longitude_reference_deg=grid.longitudes_deg)
        ax1.plot(along_km, alts_km, color=o_cmap(0.35 + 0.6 * freq_norm(solution.frequency_mhz)), linewidth=1.8, alpha=0.95)
    for solution in homed_x:
        _, _, along_km, alts_km = _ray_profile(result.scenario.tx, solution.ray, longitude_reference_deg=grid.longitudes_deg)
        ax1.plot(along_km, alts_km, color=x_cmap(0.35 + 0.6 * freq_norm(solution.frequency_mhz)), linewidth=1.7, alpha=0.92)
    ax1.scatter(
        [0.0, result.scenario.distance_km],
        [result.scenario.tx.alt_km, result.rx.alt_km],
        c=["gold", "tab:red"],
        edgecolors="black",
        linewidths=0.8,
        s=68,
        zorder=5,
    )
    ax1.set_xlim(*xlim1)
    ax1.set_ylim(*ylim1)
    ax1.set_xlabel("Ground-track distance from transmitter (km)")
    ax1.set_ylabel("Altitude (km)")
    ax1.set_title("Altitude Profile with Plasma-Frequency Swath")
    ax1.set_axisbelow(True)
    ax1.minorticks_on()
    ax1.grid(True, which="major", alpha=0.28, linewidth=0.7)
    ax1.grid(True, which="minor", alpha=0.12, linewidth=0.5, linestyle=":")
    ax1.legend(
        handles=[
            Line2D([0.0], [0.0], color=o_cmap(0.8), linewidth=2.0, label="O-mode homed"),
            Line2D([0.0], [0.0], color=x_cmap(0.8), linewidth=1.9, label="X-mode homed"),
        ],
        loc="best",
        framealpha=0.95,
    )
    cbar1 = fig.colorbar(surface1, ax=ax1, pad=0.02, label="Plasma frequency (MHz)")
    cbar1.ax.tick_params(labelsize=10)

    ax2.axhline(
        result.line_of_sight_slant_range_km,
        color="0.45",
        linestyle="--",
        linewidth=1.1,
        label="LOS slant range",
    )
    if iono_o_freqs.size:
        ax2.scatter(iono_o_freqs, iono_o_ranges, color=o_cmap(0.78), s=26, label="O-mode homed")
    if iono_x_freqs.size:
        ax2.scatter(iono_x_freqs, iono_x_ranges, color=x_cmap(0.82), s=26, label="X-mode homed")
    ax2.set_xlim(result.scenario.frequency_start_mhz, result.scenario.frequency_stop_mhz)
    ax2.set_ylim(iono_ymin - iono_margin, iono_ymax + iono_margin)
    ax2.set_xlabel("Frequency (MHz)")
    ax2.set_ylabel("Virtual range (km)")
    ax2.set_title("Sweep Ionogram")
    ax2.invert_yaxis()
    ax2.set_axisbelow(True)
    ax2.minorticks_on()
    ax2.grid(True, which="major", alpha=0.28, linewidth=0.7)
    ax2.grid(True, which="minor", alpha=0.12, linewidth=0.5, linestyle=":")
    ax2.legend(loc="best", framealpha=0.95)

    fig.suptitle(
        (
            f"Space-to-Space Example: {result.scenario.frequency_start_mhz:.1f}-{result.scenario.frequency_stop_mhz:.1f} MHz "
            f"({result.scenario.frequency_step_khz:.0f} kHz step), "
            f"{result.scenario.tx.alt_km:.0f} km to {result.rx.alt_km:.0f} km"
        ),
        fontsize=15,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.01,
        (
            f"LOS elev={result.line_of_sight_elevation_deg:.2f} deg | "
            f"O homed returns={len(homed_o)} | "
            f"X homed returns={len(homed_x)} | "
            "darker rays = higher frequency | ionogram uses homed group path"
        ),
        ha="center",
        va="bottom",
        fontsize=9,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.94))
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return path


def run_space_to_space_demo(
    scenario: SpaceToSpaceScenario | None = None,
    *,
    grid_cache_path: Path | None = None,
    rebuild_grid: bool = False,
) -> SpaceToSpaceDemoResult:
    scenario = scenario or default_scenario()
    rx = destination_point(
        scenario.tx,
        scenario.bearing_deg,
        scenario.distance_km,
        alt_km=scenario.rx_alt_km,
    )
    los_bearing_deg, los_elevation_deg, slant_range_km = _line_of_sight_angles(scenario.tx, rx)
    grid = _load_or_build_demo_grid(
        scenario,
        rx,
        grid_cache_path=grid_cache_path,
        rebuild_grid=rebuild_grid,
    )
    tracer = PointToPointRayTracer()
    sweep = _solve_frequency_sweep(
        tracer,
        tx=scenario.tx,
        rx=rx,
        grid=grid,
        frequencies_mhz=scenario.frequencies_mhz,
        nhops=scenario.nhops,
        homing_tolerance_m=scenario.homing_tolerance_m,
        start_angles_deg=(los_elevation_deg, los_bearing_deg),
    )
    return SpaceToSpaceDemoResult(
        scenario=scenario,
        rx=rx,
        line_of_sight_bearing_deg=los_bearing_deg,
        line_of_sight_elevation_deg=los_elevation_deg,
        line_of_sight_slant_range_km=slant_range_km,
        sweep=sweep,
    )


def build_parser() -> argparse.ArgumentParser:
    scenario = default_scenario()
    parser = argparse.ArgumentParser(
        description="Space-to-space example that explicitly prepares and passes a PHaRLAP ray state vector.",
    )
    parser.add_argument("--time", default=scenario.when.isoformat(), help="UTC time in ISO-8601 format.")
    parser.add_argument("--tx", nargs=3, type=float, metavar=("LAT", "LON", "ALT_KM"),
                        default=(scenario.tx.lat_deg, scenario.tx.lon_deg, scenario.tx.alt_km))
    parser.add_argument("--distance-km", type=float, default=scenario.distance_km)
    parser.add_argument("--bearing-deg", type=float, default=scenario.bearing_deg)
    parser.add_argument("--rx-alt-km", type=float, default=scenario.rx_alt_km)
    parser.add_argument("--freq-start-mhz", type=float, default=scenario.frequency_start_mhz)
    parser.add_argument("--freq-stop-mhz", type=float, default=scenario.frequency_stop_mhz)
    parser.add_argument("--freq-step-khz", type=float, default=scenario.frequency_step_khz)
    parser.add_argument("--f107", type=float, default=scenario.f107)
    parser.add_argument("--ap-daily", type=float, default=scenario.ap_daily)
    parser.add_argument("--nhops", type=int, default=scenario.nhops)
    parser.add_argument("--homing-tolerance-m", type=float, default=scenario.homing_tolerance_m)
    parser.add_argument("--grid-cache", type=Path, default=None, help="Optional saved ionosphere model cache path.")
    parser.add_argument("--rebuild-grid", action="store_true", help="Ignore any saved grid cache and rebuild the model.")
    parser.add_argument("--plot-out", type=Path, default=None, help="Optional PNG overview output path.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    scenario = SpaceToSpaceScenario(
        when=_parse_time(args.time),
        tx=GeoPoint(*args.tx),
        distance_km=float(args.distance_km),
        bearing_deg=float(args.bearing_deg),
        rx_alt_km=float(args.rx_alt_km),
        frequency_start_mhz=float(args.freq_start_mhz),
        frequency_stop_mhz=float(args.freq_stop_mhz),
        frequency_step_khz=float(args.freq_step_khz),
        f107=float(args.f107),
        ap_daily=float(args.ap_daily),
        nhops=int(args.nhops),
        homing_tolerance_m=float(args.homing_tolerance_m),
    )
    try:
        result = run_space_to_space_demo(
            scenario,
            grid_cache_path=args.grid_cache,
            rebuild_grid=bool(args.rebuild_grid),
        )
    except PyLapImportError as exc:
        raise SystemExit(str(exc)) from exc
    if args.plot_out is not None:
        plot_space_to_space_overview(
            args.plot_out.expanduser(),
            result,
            grid_cache_path=args.grid_cache,
            rebuild_grid=bool(args.rebuild_grid),
        )
    print(json.dumps(result.to_dict(), indent=2))


if __name__ == "__main__":
    main()
