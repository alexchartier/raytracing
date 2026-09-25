from __future__ import annotations

import argparse
import gc
import datetime as dt
import hashlib
import json
import math
import multiprocessing as mp
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable, Sequence

import netCDF4
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import gaussian_filter
from scipy.optimize import differential_evolution, minimize

if __package__:
    from .absorption import effective_collision_frequency
    from .geometry import GeoPoint, coerce_longitude_for_grid, great_circle_waypoints, initial_bearing_deg, llh_to_ecef, wrap_longitude
    from .grid import (
        IonosphereGrid,
        build_pyiri_grid_from_axes,
        extract_ionosphere_subgrid,
        load_ionosphere_grid_netcdf,
        save_ionosphere_grid_netcdf,
    )
    from .tracer import PointToPointRayTracer, PyLapImportError
else:
    from absorption import effective_collision_frequency
    from geometry import GeoPoint, coerce_longitude_for_grid, great_circle_waypoints, initial_bearing_deg, llh_to_ecef, wrap_longitude
    from grid import (
        IonosphereGrid,
        build_pyiri_grid_from_axes,
        extract_ionosphere_subgrid,
        load_ionosphere_grid_netcdf,
        save_ionosphere_grid_netcdf,
    )
    from tracer import PointToPointRayTracer, PyLapImportError


C_M_PER_S = 299_792_458.0
DEFAULT_AMPERE_FILE = Path("~/Downloads/ampere.20100101.rawdB.nc").expanduser()


@dataclass(frozen=True)
class IonosphereFitParams:
    density_scale: float = 1.0
    hmf2_shift_km: float = 0.0
    wave_amplitude_fraction: float = 0.10
    wave_phase_rad: float = 0.0
    wave_bearing_deg: float = 35.0

    def to_vector(self, parameter_names: Sequence[str] | None = None) -> np.ndarray:
        names = _fit_parameter_names() if parameter_names is None else tuple(parameter_names)
        return np.array([getattr(self, name) for name in names], dtype=float)

    @classmethod
    def from_vector(
        cls,
        values: Sequence[float],
        parameter_names: Sequence[str] | None = None,
        *,
        base: IonosphereFitParams | None = None,
    ) -> IonosphereFitParams:
        names = _fit_parameter_names() if parameter_names is None else tuple(parameter_names)
        params = {
            "density_scale": 1.0,
            "hmf2_shift_km": 0.0,
            "wave_amplitude_fraction": 0.10,
            "wave_phase_rad": 0.0,
            "wave_bearing_deg": 35.0,
        }
        if base is not None:
            params.update(vars(base))
        for name, value in zip(names, values):
            params[name] = float(value)
        return cls(**params)


@dataclass(frozen=True)
class TopsideInverseConfig:
    ampere_file: Path = DEFAULT_AMPERE_FILE
    center_utc: dt.datetime = dt.datetime(2010, 1, 1, 12, 0, 0)
    planes: tuple[int, ...] = (0, 1, 3, 5)
    spoof_altitude_km: float = 800.0
    oblique_separation_km: float = 600.0
    epoch_offsets: tuple[int, ...] = (-1, 0, 1)
    frequencies_mhz: tuple[float, ...] = (6.0, 8.0, 10.0, 12.0, 14.0)
    range_min_km: float = 150.0
    range_max_km: float = 2600.0
    range_bin_km: float = 25.0
    blur_sigma_frequency_bins: float = 0.35
    blur_sigma_range_bins: float = 0.75
    homing_tolerance_m: float = 1_000.0
    seed_max_candidates_per_frequency: int = 8
    homed_max_returns_per_frequency: int = 10
    vertical_elevation_min_deg: float = -88.0
    vertical_elevation_max_deg: float = -12.0
    vertical_elevation_count: int = 24
    vertical_azimuth_step_deg: float = 30.0
    oblique_elevation_count: int = 31
    oblique_bearing_count: int = 15
    grid_lat_step_deg: float = 0.5
    grid_lon_step_deg: float = 1.0
    grid_lat_margin_deg: float = 5.0
    grid_lon_margin_deg: float = 5.0
    vertical_subgrid_margin_deg: float = 24.0
    oblique_subgrid_margin_deg: float = 8.0
    grid_alt_min_km: float = 60.0
    grid_alt_max_km: float = 900.0
    grid_alt_step_km: float = 5.0
    f107: float = 120.0
    ap_daily: float = 8.0
    d_region_model: str = "fpt2018"
    wave_horizontal_wavelength_km: float = 450.0
    wave_vertical_center_km: float = 280.0
    wave_vertical_sigma_km: float = 60.0
    initial_params: IonosphereFitParams = IonosphereFitParams()
    fit_parameter_names: tuple[str, ...] = (
        "density_scale",
        "hmf2_shift_km",
        "wave_amplitude_fraction",
        "wave_phase_rad",
        "wave_bearing_deg",
    )
    solver_maxiter: int = 10
    solver_popsize: int = 8
    solver_seed: int = 0
    grid_cache_path: Path | None = None
    rebuild_grid: bool = False
    truth_params: IonosphereFitParams = IonosphereFitParams(
        density_scale=1.035,
        hmf2_shift_km=15.0,
        wave_amplitude_fraction=0.12,
        wave_phase_rad=0.55,
        wave_bearing_deg=32.0,
    )


@dataclass(frozen=True)
class SpoofSatelliteTrack:
    plane: int
    pseudo_sv: int
    times_utc: tuple[dt.datetime, ...]
    positions: tuple[GeoPoint, ...]


@dataclass(frozen=True)
class TopsideCase:
    name: str
    kind: str
    plane: int
    pseudo_sv: int
    times_utc: tuple[dt.datetime, ...]
    tx_points: tuple[GeoPoint, ...]
    rx_points: tuple[GeoPoint, ...]
    fan_elevations_deg: np.ndarray
    fan_bearings_deg: np.ndarray


@dataclass(frozen=True)
class CaseObservables:
    name: str
    kind: str
    total_image: np.ndarray
    o_image: np.ndarray
    x_image: np.ndarray
    ox_split_image: np.ndarray
    doppler_image_hz: np.ndarray
    doppler_weight: np.ndarray


@dataclass(frozen=True)
class HomedRayReturn:
    frequency_mhz: float
    ox_mode: int
    ray: object
    miss_m: float
    group_range_km: float
    absorption_db: float
    doppler_hz: float


@dataclass(frozen=True)
class RayProfileCurve:
    frequency_mhz: float
    ox_mode: int
    miss_m: float
    group_range_km: float
    along_km: np.ndarray
    alt_km: np.ndarray


@dataclass(frozen=True)
class FrequencyPlotSlice:
    frequency_mhz: float
    o_image: np.ndarray
    x_image: np.ndarray
    doppler_num: np.ndarray
    doppler_den: np.ndarray
    o_curves: tuple[RayProfileCurve, ...]
    x_curves: tuple[RayProfileCurve, ...]


@dataclass(frozen=True)
class SyntheticDataset:
    frequencies_mhz: np.ndarray
    range_edges_km: np.ndarray
    range_centers_km: np.ndarray
    cases: tuple[CaseObservables, ...]


@dataclass(frozen=True)
class InverseProblem:
    config: TopsideInverseConfig
    cases: tuple[TopsideCase, ...]
    global_background_grid: IonosphereGrid | None
    background_grids: tuple[IonosphereGrid, ...]
    wave_origin: GeoPoint
    frequencies_mhz: np.ndarray
    range_edges_km: np.ndarray
    range_centers_km: np.ndarray


@dataclass(frozen=True)
class FitResult:
    truth_params: IonosphereFitParams
    fitted_params: IonosphereFitParams
    truth_cost: float
    fitted_cost: float
    solver_success: bool
    solver_message: str
    iterations: int
    evaluations: int

    def to_dict(self) -> dict:
        return {
            "truth_params": vars(self.truth_params),
            "fitted_params": vars(self.fitted_params),
            "truth_cost": self.truth_cost,
            "fitted_cost": self.fitted_cost,
            "solver_success": self.solver_success,
            "solver_message": self.solver_message,
            "iterations": self.iterations,
            "evaluations": self.evaluations,
        }


@dataclass(frozen=True)
class InverseDemoResult:
    problem: InverseProblem
    observed: SyntheticDataset
    fitted: SyntheticDataset
    fit: FitResult

    def to_dict(self) -> dict:
        return {
            "config": {
                "center_utc": self.problem.config.center_utc.isoformat(),
                "planes": list(self.problem.config.planes),
                "spoof_altitude_km": self.problem.config.spoof_altitude_km,
                "oblique_separation_km": self.problem.config.oblique_separation_km,
                "frequencies_mhz": list(self.problem.frequencies_mhz),
            },
            "cases": [
                {
                    "name": case.name,
                    "kind": case.kind,
                    "plane": case.plane,
                    "pseudo_sv": case.pseudo_sv,
                    "times_utc": [when.isoformat() for when in case.times_utc],
                    "tx": [vars(point) for point in case.tx_points],
                    "rx": [vars(point) for point in case.rx_points],
                }
                for case in self.problem.cases
            ],
            "fit": self.fit.to_dict(),
        }


def _subset_problem_cases(problem: InverseProblem, case_names: Sequence[str] | None = None) -> InverseProblem:
    if case_names is None:
        return problem
    requested = tuple(case_names)
    if not requested:
        return problem
    selected_indices = [index for index, case in enumerate(problem.cases) if case.name in requested]
    if not selected_indices:
        raise ValueError(f"unknown case names: {', '.join(requested)}")
    missing = [name for name in requested if all(case.name != name for case in problem.cases)]
    if missing:
        raise ValueError(f"unknown case names: {', '.join(missing)}")
    return replace(
        problem,
        cases=tuple(problem.cases[index] for index in selected_indices),
        background_grids=tuple(problem.background_grids[index] for index in selected_indices),
    )


def _empty_case_observables(problem: InverseProblem, case: TopsideCase) -> CaseObservables:
    shape = (problem.frequencies_mhz.size, problem.range_centers_km.size)
    zeros = np.zeros(shape, dtype=float)
    return CaseObservables(
        name=case.name,
        kind=case.kind,
        total_image=zeros.copy(),
        o_image=zeros.copy(),
        x_image=zeros.copy(),
        ox_split_image=zeros.copy(),
        doppler_image_hz=zeros.copy(),
        doppler_weight=zeros.copy(),
    )


def _load_fit_result_json(path: Path) -> FitResult:
    payload = json.loads(path.read_text())
    fit_payload = payload.get("fit", payload)
    return FitResult(
        truth_params=IonosphereFitParams(**fit_payload["truth_params"]),
        fitted_params=IonosphereFitParams(**fit_payload["fitted_params"]),
        truth_cost=float(fit_payload.get("truth_cost", float("nan"))),
        fitted_cost=float(fit_payload.get("fitted_cost", float("nan"))),
        solver_success=bool(fit_payload.get("solver_success", True)),
        solver_message=str(fit_payload.get("solver_message", "")),
        iterations=int(fit_payload.get("iterations", 0)),
        evaluations=int(fit_payload.get("evaluations", 0)),
    )


def _fit_parameter_names() -> tuple[str, ...]:
    return (
        "density_scale",
        "hmf2_shift_km",
        "wave_amplitude_fraction",
        "wave_phase_rad",
        "wave_bearing_deg",
    )


def _fit_parameter_bounds() -> dict[str, tuple[float, float]]:
    return {
        "density_scale": (0.85, 1.15),
        "hmf2_shift_km": (-40.0, 40.0),
        "wave_amplitude_fraction": (0.0, 0.25),
        "wave_phase_rad": (-math.pi, math.pi),
        "wave_bearing_deg": (-90.0, 90.0),
    }


def _parse_time(text: str) -> dt.datetime:
    when = dt.datetime.fromisoformat(text)
    if when.tzinfo is not None:
        return when.astimezone(dt.timezone.utc).replace(tzinfo=None)
    return when


def _gmst_from_jd(jd_ut1: float) -> float:
    t = (jd_ut1 - 2451545.0) / 36525.0
    gmst_deg = (
        280.46061837
        + 360.98564736629 * (jd_ut1 - 2451545.0)
        + 0.000387933 * t * t
        - t * t * t / 38710000.0
    )
    return math.radians(gmst_deg % 360.0)


def _teme_to_ecef(r_km: np.ndarray, jd_ut1: float) -> np.ndarray:
    theta = _gmst_from_jd(jd_ut1)
    c = math.cos(theta)
    s = math.sin(theta)
    r_m = np.asarray(r_km, dtype=np.float64) * 1000.0
    return np.array(
        [
            c * r_m[0] + s * r_m[1],
            -s * r_m[0] + c * r_m[1],
            r_m[2],
        ],
        dtype=np.float64,
    )


def _ecef_to_geo(point_xyz_m: np.ndarray) -> GeoPoint:
    x, y, z = (float(value) for value in np.asarray(point_xyz_m, dtype=np.float64))
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    wgs84_a_m = 6378137.0
    wgs84_f = 1.0 / 298.257223563
    wgs84_e2 = wgs84_f * (2.0 - wgs84_f)
    lat = math.atan2(z, p * (1.0 - wgs84_e2))
    alt_m = 0.0
    for _ in range(7):
        sin_lat = math.sin(lat)
        cos_lat = math.cos(lat)
        prime_vertical = wgs84_a_m / math.sqrt(1.0 - wgs84_e2 * sin_lat * sin_lat)
        safe_cos = cos_lat if abs(cos_lat) > 1e-12 else math.copysign(1e-12, cos_lat if cos_lat != 0.0 else 1.0)
        alt_m = p / safe_cos - prime_vertical
        denom = max(prime_vertical + alt_m, 1.0)
        lat = math.atan2(z, p * (1.0 - wgs84_e2 * prime_vertical / denom))
    return GeoPoint(math.degrees(lat), math.degrees(lon), alt_m / 1000.0)


def _to_utc(year: int, doy: int, ut_hours: float) -> dt.datetime:
    start = dt.datetime(int(year), 1, 1) + dt.timedelta(days=int(doy) - 1)
    return start + dt.timedelta(hours=float(ut_hours))


def _surface_distance_km(start: GeoPoint, end: GeoPoint) -> float:
    lat1 = math.radians(start.lat_deg)
    lat2 = math.radians(end.lat_deg)
    dlat = lat2 - lat1
    dlon = math.radians(((end.lon_deg - start.lon_deg + 180.0) % 360.0) - 180.0)
    a = math.sin(dlat / 2.0) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2
    return 6371.0088 * 2.0 * math.asin(min(1.0, math.sqrt(max(a, 0.0))))


def _coerce_longitudes_near_reference(longitudes_deg: np.ndarray, reference_deg: float) -> np.ndarray:
    longitudes = np.asarray(longitudes_deg, dtype=np.float64)
    candidates = np.stack((longitudes - 360.0, longitudes, longitudes + 360.0), axis=1)
    index = np.argmin(np.abs(candidates - float(reference_deg)), axis=1)
    return candidates[np.arange(longitudes.size), index]


def _local_east_north_km(origin: GeoPoint, lat_deg: np.ndarray | float, lon_deg: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
    lat = np.asarray(lat_deg, dtype=np.float64)
    lon = np.asarray(lon_deg, dtype=np.float64)
    origin_lat_rad = math.radians(origin.lat_deg)
    dlat_rad = np.deg2rad(lat - origin.lat_deg)
    flat_lon = np.ravel(lon)
    dlon_rad = np.deg2rad(
        np.asarray([wrap_longitude(float(value) - origin.lon_deg) for value in flat_lon], dtype=np.float64)
    ).reshape(lon.shape)
    north_km = 6371.0088 * dlat_rad
    east_km = 6371.0088 * math.cos(origin_lat_rad) * dlon_rad
    return east_km, north_km


def _along_track_coordinate_km(origin: GeoPoint, bearing_deg: float, lat_deg: np.ndarray | float, lon_deg: np.ndarray | float) -> np.ndarray:
    east_km, north_km = _local_east_north_km(origin, lat_deg, lon_deg)
    az = math.radians(bearing_deg)
    return east_km * math.sin(az) + north_km * math.cos(az)


def _centers_to_edges(values: np.ndarray) -> np.ndarray:
    centers = np.asarray(values, dtype=float)
    if centers.size == 1:
        return np.array([centers[0] - 0.5, centers[0] + 0.5], dtype=float)
    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def make_frequency_axis_mhz(start_mhz: float, stop_mhz: float, step_khz: float) -> np.ndarray:
    if stop_mhz < start_mhz:
        raise ValueError("stop_mhz must be greater than or equal to start_mhz")
    step_mhz = max(float(step_khz) / 1000.0, 1e-6)
    count = int(math.floor((float(stop_mhz) - float(start_mhz)) / step_mhz + 1e-9)) + 1
    axis = float(start_mhz) + step_mhz * np.arange(max(count, 1), dtype=float)
    return axis[axis <= float(stop_mhz) + 1e-9]


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


def _return_leg_distance(ray_path: dict, tx: GeoPoint, rx: GeoPoint) -> tuple[float, float | None, float | None]:
    heights = np.asarray(ray_path.get("height", []), dtype=float)
    lats = np.asarray(ray_path.get("lat", []), dtype=float)
    lons = np.asarray(ray_path.get("lon", []), dtype=float)
    valid = np.isfinite(heights) & np.isfinite(lats) & np.isfinite(lons) & (heights < 1e40)
    if np.count_nonzero(valid) < 2:
        return math.inf, None, None

    heights = heights[valid]
    lats = lats[valid]
    lons = lons[valid]
    ray_xyz = llh_to_ecef(lats, lons, heights * 1000.0)
    point_xyz = np.asarray(llh_to_ecef(rx.lat_deg, rx.lon_deg, rx.alt_km * 1000.0), dtype=float)
    if point_xyz.ndim > 1:
        point_xyz = point_xyz[0]

    group_range = np.asarray(ray_path.get("group_range", np.full_like(heights, np.nan)), dtype=float)[valid]
    geom_distance = np.asarray(ray_path.get("geometric_distance", np.full_like(heights, np.nan)), dtype=float)[valid]
    absorption = np.asarray(ray_path.get("absorption", np.full_like(heights, np.nan)), dtype=float)[valid]

    start_index = 1
    if tx.alt_km > 400.0 and rx.alt_km > 400.0:
        threshold_alt_km = min(tx.alt_km, rx.alt_km) - 20.0
        below = np.flatnonzero(heights < threshold_alt_km)
        if below.size > 0:
            start_index = int(below[0])
    if start_index >= len(heights) - 1:
        return math.inf, None, None

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
        return math.inf, None, None
    best_index = int(np.argmin(distances))
    best_distance = float(distances[best_index])
    segment_index = start_index + best_index
    best_group_km: float | None = None
    best_absorption_db: float | None = None
    if np.isfinite(group_range[segment_index]) and np.isfinite(group_range[segment_index + 1]):
        best_group_km = float(group_range[segment_index] + t[best_index] * (group_range[segment_index + 1] - group_range[segment_index]))
    elif np.isfinite(geom_distance[segment_index]) and np.isfinite(geom_distance[segment_index + 1]):
        best_group_km = float(geom_distance[segment_index] + t[best_index] * (geom_distance[segment_index + 1] - geom_distance[segment_index]))
    if np.isfinite(absorption[segment_index]) and np.isfinite(absorption[segment_index + 1]):
        best_absorption_db = float(absorption[segment_index] + t[best_index] * (absorption[segment_index + 1] - absorption[segment_index]))
    return best_distance, best_group_km, best_absorption_db


def _doppler_from_summary(ray) -> float:
    value = ray.summary.get("Doppler_shift", np.nan)
    try:
        value = float(np.asarray(value).ravel()[0])
    except Exception:
        return float("nan")
    return value if math.isfinite(value) else float("nan")


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


def _topside_search_arrays(tx: GeoPoint, rx: GeoPoint, config: TopsideInverseConfig) -> tuple[np.ndarray, np.ndarray]:
    gc_bearing_deg, gc_elevation_deg, slant_range_km = _line_of_sight_angles(tx, rx)
    midpoint = GeoPoint(
        0.5 * (tx.lat_deg + rx.lat_deg),
        0.5 * (tx.lon_deg + rx.lon_deg),
        0.0,
    )
    tangent_elevation_deg = _line_of_sight_angles(tx, midpoint)[1]
    if not math.isfinite(tangent_elevation_deg):
        tangent_elevation_deg = gc_elevation_deg - 40.0
    elevation_min = max(-89.0, min(tangent_elevation_deg, gc_elevation_deg))
    elevation_max = min(89.0, max(tangent_elevation_deg, gc_elevation_deg))
    elevations_deg = np.linspace(elevation_min, elevation_max, max(config.oblique_elevation_count, 3), dtype=float)
    az_half_width_deg = max(3.0, float(math.ceil(30.0 * math.sqrt(100.0 / max(slant_range_km, 100.0)))))
    bearings_deg = np.linspace(
        gc_bearing_deg - az_half_width_deg,
        gc_bearing_deg + az_half_width_deg,
        max(config.oblique_bearing_count, 3),
        dtype=float,
    )
    return elevations_deg, bearings_deg


def _vertical_search_arrays(config: TopsideInverseConfig) -> tuple[np.ndarray, np.ndarray]:
    elevations_deg = np.linspace(
        config.vertical_elevation_min_deg,
        config.vertical_elevation_max_deg,
        max(config.vertical_elevation_count, 3),
        dtype=float,
    )
    azimuths_deg = np.arange(0.0, 360.0, max(config.vertical_azimuth_step_deg, 1.0), dtype=float)
    if azimuths_deg.size < 4:
        azimuths_deg = np.linspace(0.0, 330.0, 12, dtype=float)
    return elevations_deg, azimuths_deg


def _fan_mesh(elevations_deg: np.ndarray, bearings_deg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    elevation_grid, bearing_grid = np.meshgrid(elevations_deg, bearings_deg, indexing="xy")
    return elevation_grid.ravel(), bearing_grid.ravel()


def load_ampere_tracks(config: TopsideInverseConfig) -> tuple[SpoofSatelliteTrack, ...]:
    ampere_file = config.ampere_file.expanduser()
    with netCDF4.Dataset(ampere_file) as dataset:
        time_hours = np.asarray(dataset.variables["time"][:], dtype=float)
        planes = np.asarray(dataset.variables["plane_num"][:], dtype=int)
        pseudo_sv = np.asarray(dataset.variables["pseudo_sv_num"][:], dtype=int)
        pos_eci_m = np.asarray(dataset.variables["pos_eci"][:], dtype=float)
        years = np.asarray(dataset.variables["year"][:], dtype=int)
        doys = np.asarray(dataset.variables["doy"][:], dtype=int)

    center_hour = (
        config.center_utc.hour
        + config.center_utc.minute / 60.0
        + config.center_utc.second / 3600.0
        + config.center_utc.microsecond / 3.6e9
    )
    tracks: list[SpoofSatelliteTrack] = []
    for plane in config.planes:
        plane_indices = np.flatnonzero(planes == int(plane))
        if plane_indices.size == 0:
            raise RuntimeError(f"AMPERE file does not contain plane {plane}")
        center_index = int(plane_indices[np.argmin(np.abs(time_hours[plane_indices] - center_hour))])
        sv_id = int(pseudo_sv[center_index])
        track_indices = np.flatnonzero((planes == int(plane)) & (pseudo_sv == sv_id))
        order = np.argsort(time_hours[track_indices])
        track_indices = track_indices[order]
        track_times: list[dt.datetime] = []
        track_positions: list[GeoPoint] = []
        for raw_index in track_indices:
            when_utc = _to_utc(int(years[raw_index]), int(doys[raw_index]), float(time_hours[raw_index]))
            jd = when_utc.timestamp() / 86400.0 + 2440587.5
            ecef_m = _teme_to_ecef(pos_eci_m[raw_index] / 1000.0, jd)
            geo = _ecef_to_geo(ecef_m)
            track_times.append(when_utc)
            track_positions.append(GeoPoint(geo.lat_deg, geo.lon_deg, config.spoof_altitude_km))
        tracks.append(
            SpoofSatelliteTrack(
                plane=int(plane),
                pseudo_sv=sv_id,
                times_utc=tuple(track_times),
                positions=tuple(track_positions),
            )
        )
    return tuple(tracks)


def _find_oblique_offset(track: SpoofSatelliteTrack, center_index: int, desired_distance_km: float) -> int:
    tx = track.positions[center_index]
    distances = np.asarray(
        [
            _surface_distance_km(GeoPoint(tx.lat_deg, tx.lon_deg, 0.0), GeoPoint(point.lat_deg, point.lon_deg, 0.0))
            for point in track.positions[center_index:]
        ],
        dtype=float,
    )
    positive = np.flatnonzero(distances > 0.0)
    if positive.size == 0:
        raise RuntimeError("track does not contain a forward oblique separation")
    local_index = int(positive[np.argmin(np.abs(distances[positive] - desired_distance_km))])
    return local_index


def build_topside_cases(config: TopsideInverseConfig) -> tuple[TopsideCase, ...]:
    tracks = load_ampere_tracks(config)
    offsets = np.asarray(config.epoch_offsets, dtype=int)
    if offsets.size != 3 or tuple(offsets) != (-1, 0, 1):
        raise ValueError("epoch_offsets must currently be exactly (-1, 0, 1) for centered Doppler differencing")

    center_time = config.center_utc
    cases: list[TopsideCase] = []
    for track in tracks:
        track_times = np.asarray(track.times_utc, dtype=object)
        center_index = int(
            np.argmin(np.abs(np.asarray([(when - center_time).total_seconds() for when in track_times], dtype=float)))
        )
        oblique_offset = _find_oblique_offset(track, center_index, config.oblique_separation_km)
        tx_indices = center_index + offsets
        rx_indices_oblique = tx_indices + oblique_offset
        if int(np.min(tx_indices)) < 0 or int(np.max(rx_indices_oblique)) >= len(track.positions):
            raise RuntimeError(f"track for plane {track.plane} does not have enough samples around {center_time.isoformat()}")

        tx_points = tuple(track.positions[int(index)] for index in tx_indices)
        tx_times = tuple(track.times_utc[int(index)] for index in tx_indices)
        vertical_elevs, vertical_bears = _vertical_search_arrays(config)
        oblique_elevs, oblique_bears = _topside_search_arrays(tx_points[1], track.positions[int(rx_indices_oblique[1])], config)
        vertical_fan = _fan_mesh(vertical_elevs, vertical_bears)
        oblique_fan = _fan_mesh(oblique_elevs, oblique_bears)

        cases.append(
            TopsideCase(
                name=f"plane{track.plane:02d}_sv{track.pseudo_sv:03d}_vertical",
                kind="vertical",
                plane=track.plane,
                pseudo_sv=track.pseudo_sv,
                times_utc=tx_times,
                tx_points=tx_points,
                rx_points=tx_points,
                fan_elevations_deg=vertical_fan[0],
                fan_bearings_deg=vertical_fan[1],
            )
        )
        cases.append(
            TopsideCase(
                name=f"plane{track.plane:02d}_sv{track.pseudo_sv:03d}_oblique",
                kind="oblique",
                plane=track.plane,
                pseudo_sv=track.pseudo_sv,
                times_utc=tx_times,
                tx_points=tx_points,
                rx_points=tuple(track.positions[int(index)] for index in rx_indices_oblique),
                fan_elevations_deg=oblique_fan[0],
                fan_bearings_deg=oblique_fan[1],
            )
        )
    return tuple(cases)


def _global_grid_cache_key(config: TopsideInverseConfig, cases: Sequence[TopsideCase], lat_axis: np.ndarray, lon_axis: np.ndarray) -> str:
    payload = {
        "center_utc": config.center_utc.isoformat(),
        "cases": [case.name for case in cases],
        "lat_axis": [float(lat_axis[0]), float(lat_axis[-1]), int(lat_axis.size)],
        "lon_axis": [float(lon_axis[0]), float(lon_axis[-1]), int(lon_axis.size)],
        "grid": {
            "lat_step": config.grid_lat_step_deg,
            "lon_step": config.grid_lon_step_deg,
            "alt_min": config.grid_alt_min_km,
            "alt_max": config.grid_alt_max_km,
            "alt_step": config.grid_alt_step_km,
            "f107": config.f107,
            "ap_daily": config.ap_daily,
            "d_region_model": config.d_region_model,
        },
    }
    return hashlib.sha1(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]


def _case_support_points(case: TopsideCase) -> tuple[np.ndarray, np.ndarray]:
    latitudes: list[float] = []
    longitudes: list[float] = []
    wrapped = np.asarray(
        [wrap_longitude(point.lon_deg) for point in (*case.tx_points, *case.rx_points)],
        dtype=float,
    )
    reference_lon = float(np.mean(wrapped))
    for tx, rx in zip(case.tx_points, case.rx_points):
        tx_lon = _coerce_longitudes_near_reference(np.array([wrap_longitude(tx.lon_deg)], dtype=float), reference_lon)[0]
        rx_lon = _coerce_longitudes_near_reference(np.array([wrap_longitude(rx.lon_deg)], dtype=float), reference_lon)[0]
        latitudes.extend((tx.lat_deg, rx.lat_deg))
        longitudes.extend((tx_lon, rx_lon))
        if case.kind == "oblique":
            path_lats, path_lons = great_circle_waypoints(
                GeoPoint(tx.lat_deg, tx.lon_deg, 0.0),
                GeoPoint(rx.lat_deg, rx.lon_deg, 0.0),
                count=64,
            )
            path_lons = _coerce_longitudes_near_reference(np.asarray(path_lons, dtype=float), reference_lon)
            latitudes.extend(float(value) for value in path_lats)
            longitudes.extend(float(value) for value in path_lons)
    return np.asarray(latitudes, dtype=float), np.asarray(longitudes, dtype=float)


def _case_grid_axes(case: TopsideCase, config: TopsideInverseConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    latitudes, longitudes = _case_support_points(case)
    margin_deg = config.vertical_subgrid_margin_deg if case.kind == "vertical" else config.oblique_subgrid_margin_deg
    lat_step = max(config.grid_lat_step_deg, 1e-6)
    lon_step = max(config.grid_lon_step_deg, 1e-6)
    alt_step = max(config.grid_alt_step_km, 1e-6)
    lat_axis = np.array([], dtype=float)
    lon_axis = np.array([], dtype=float)
    current_margin_deg = margin_deg
    while True:
        lat_min = max(-90.0, float(np.min(latitudes)) - current_margin_deg)
        lat_max = min(90.0, float(np.max(latitudes)) + current_margin_deg)
        lon_min = float(np.min(longitudes)) - current_margin_deg
        lon_max = float(np.max(longitudes)) + current_margin_deg
        lat_axis = np.arange(
            math.floor(lat_min / lat_step) * lat_step,
            math.ceil(lat_max / lat_step) * lat_step + 0.5 * lat_step,
            lat_step,
            dtype=float,
        )
        lon_start = math.floor(lon_min / lon_step) * lon_step
        lon_stop = math.ceil(lon_max / lon_step) * lon_step
        while lon_start < -180.0:
            lon_start += 360.0
            lon_stop += 360.0
        while lon_start > 180.0:
            lon_start -= 360.0
            lon_stop -= 360.0
        lon_axis = np.arange(
            lon_start,
            lon_stop + 0.5 * lon_step,
            lon_step,
            dtype=float,
        )
        if lat_axis.size <= 101 and lon_axis.size <= 101:
            break
        if case.kind != "vertical" or current_margin_deg <= max(lat_step, lon_step):
            break
        current_margin_deg = max(current_margin_deg - max(lat_step, lon_step), 0.0)
    alt_axis = np.arange(
        config.grid_alt_min_km,
        config.grid_alt_max_km + 0.5 * alt_step,
        alt_step,
        dtype=float,
    )
    if lat_axis.size > 101 or lon_axis.size > 101 or alt_axis.size > 201:
        raise RuntimeError(
            f"{case.name} envelope exceeds PHaRLAP grid limits at requested resolution: "
            f"{lat_axis.size} lat x {lon_axis.size} lon x {alt_axis.size} alt"
        )
    return lat_axis, lon_axis, alt_axis


def _global_grid_axes(cases: Sequence[TopsideCase], config: TopsideInverseConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    support_lats = []
    for case in cases:
        latitudes, _ = _case_support_points(case)
        support_lats.extend(float(value) for value in latitudes)
    if not support_lats:
        raise RuntimeError("no topside cases were built")
    lat_step = max(config.grid_lat_step_deg, 1e-6)
    lon_step = max(config.grid_lon_step_deg, 1e-6)
    alt_step = max(config.grid_alt_step_km, 1e-6)
    lat_pad = max(config.grid_lat_margin_deg, config.vertical_subgrid_margin_deg, config.oblique_subgrid_margin_deg)
    lat_min = max(-90.0, min(support_lats) - lat_pad)
    lat_max = min(90.0, max(support_lats) + lat_pad)
    lat_axis = np.arange(
        math.floor(lat_min / lat_step) * lat_step,
        math.ceil(lat_max / lat_step) * lat_step + 0.5 * lat_step,
        lat_step,
        dtype=float,
    )
    lon_axis = np.arange(-180.0, 180.0 + 0.5 * lon_step, lon_step, dtype=float)
    alt_axis = np.arange(
        config.grid_alt_min_km,
        config.grid_alt_max_km + 0.5 * alt_step,
        alt_step,
        dtype=float,
    )
    if lat_axis.size > 701 or lon_axis.size > 701 or alt_axis.size > 401:
        raise RuntimeError(
            "requested global background grid exceeds PHaRLAP/PyIRI limits: "
            f"{lat_axis.size} lat x {lon_axis.size} lon x {alt_axis.size} alt"
        )
    return lat_axis, lon_axis, alt_axis


def _default_global_grid_cache_path(config: TopsideInverseConfig, cases: Sequence[TopsideCase], lat_axis: np.ndarray, lon_axis: np.ndarray) -> Path:
    key = _global_grid_cache_key(config, cases, lat_axis, lon_axis)
    return Path(".cache") / "multisat_topside_inverse" / f"global_background_{key}.nc"


def build_global_background_grid(config: TopsideInverseConfig, cases: Sequence[TopsideCase]) -> IonosphereGrid:
    lat_axis, lon_axis, alt_axis = _global_grid_axes(cases, config)
    if config.grid_cache_path is None:
        cache_path = _default_global_grid_cache_path(config, cases, lat_axis, lon_axis)
    else:
        base = config.grid_cache_path.expanduser()
        if base.suffix.lower() == ".nc":
            cache_path = base
        else:
            cache_path = base / _default_global_grid_cache_path(config, cases, lat_axis, lon_axis).name
    if cache_path.exists() and not config.rebuild_grid:
        return load_ionosphere_grid_netcdf(cache_path)
    grid = build_pyiri_grid_from_axes(
        config.center_utc,
        lat_axis,
        lon_axis,
        alt_axis,
        f107=config.f107,
        ap_daily=config.ap_daily,
        d_region_model=config.d_region_model,
    )
    save_ionosphere_grid_netcdf(cache_path, grid)
    return grid


def extract_case_background_grid(global_grid: IonosphereGrid, case: TopsideCase, config: TopsideInverseConfig) -> IonosphereGrid:
    lat_axis, lon_axis, alt_axis = _case_grid_axes(case, config)
    return extract_ionosphere_subgrid(global_grid, lat_axis, lon_axis, alt_axis)


def _mean_ground_point(cases: Sequence[TopsideCase]) -> GeoPoint:
    latitudes = np.asarray([case.tx_points[1].lat_deg for case in cases], dtype=float)
    wrapped_longitudes = np.asarray([wrap_longitude(case.tx_points[1].lon_deg) for case in cases], dtype=float)
    reference_lon = float(np.mean(wrapped_longitudes))
    longitudes = _coerce_longitudes_near_reference(wrapped_longitudes, reference_lon)
    return GeoPoint(float(np.mean(latitudes)), float(np.mean(longitudes)), 0.0)


def build_inverse_problem(config: TopsideInverseConfig) -> InverseProblem:
    cases = build_topside_cases(config)
    global_background_grid = build_global_background_grid(config, cases)
    background_grids = tuple(extract_case_background_grid(global_background_grid, case, config) for case in cases)
    range_edges_km = np.arange(
        config.range_min_km,
        config.range_max_km + 0.5 * config.range_bin_km,
        config.range_bin_km,
        dtype=float,
    )
    range_centers_km = 0.5 * (range_edges_km[:-1] + range_edges_km[1:])
    return InverseProblem(
        config=config,
        cases=cases,
        global_background_grid=global_background_grid,
        background_grids=background_grids,
        wave_origin=_mean_ground_point(cases),
        frequencies_mhz=np.asarray(config.frequencies_mhz, dtype=float),
        range_edges_km=range_edges_km,
        range_centers_km=range_centers_km,
    )


def _shift_profiles_in_altitude(profiles: np.ndarray, altitudes_km: np.ndarray, shift_km: float) -> np.ndarray:
    if abs(shift_km) < 1e-9:
        return profiles.copy()
    shifted = np.empty_like(profiles)
    source_alts = np.asarray(altitudes_km, dtype=float) - float(shift_km)
    for lat_index in range(profiles.shape[1]):
        for lon_index in range(profiles.shape[2]):
            shifted[:, lat_index, lon_index] = np.interp(
                source_alts,
                altitudes_km,
                profiles[:, lat_index, lon_index],
                left=float(profiles[0, lat_index, lon_index]),
                right=float(profiles[-1, lat_index, lon_index]),
            )
    return shifted


def _apply_fit_params(problem: InverseProblem, params: IonosphereFitParams) -> IonosphereGrid:
    raise RuntimeError("_apply_fit_params now requires an explicit background grid")


def _apply_fit_params_to_grid(
    problem: InverseProblem,
    background_grid: IonosphereGrid,
    params: IonosphereFitParams,
) -> IonosphereGrid:
    base_grid = background_grid
    density_lat_lon_alt_m3 = np.asarray(base_grid.iono_en_grid, dtype=float) * 1e6
    density_alt_lat_lon_m3 = np.transpose(density_lat_lon_alt_m3, (2, 0, 1))
    shifted_density_alt_lat_lon_m3 = _shift_profiles_in_altitude(
        density_alt_lat_lon_m3,
        np.asarray(base_grid.altitudes_km, dtype=float),
        params.hmf2_shift_km,
    )

    lon_mesh, lat_mesh = np.meshgrid(base_grid.longitudes_deg, base_grid.latitudes_deg, indexing="xy")
    phase_coordinate_km = _along_track_coordinate_km(problem.wave_origin, params.wave_bearing_deg, lat_mesh, lon_mesh)
    vertical = np.exp(
        -0.5 * ((np.asarray(base_grid.altitudes_km, dtype=float) - problem.config.wave_vertical_center_km) / max(problem.config.wave_vertical_sigma_km, 1e-6)) ** 2
    )
    phase = (
        2.0
        * math.pi
        * phase_coordinate_km[None, :, :]
        / max(problem.config.wave_horizontal_wavelength_km, 1e-6)
        + params.wave_phase_rad
    )
    perturbation = params.wave_amplitude_fraction * vertical[:, None, None] * np.cos(phase)
    perturbed_density_alt_lat_lon_m3 = np.maximum(
        shifted_density_alt_lat_lon_m3 * max(params.density_scale, 1e-6) * (1.0 + perturbation),
        1e-6,
    )
    perturbed_density_lat_lon_alt_m3 = np.transpose(perturbed_density_alt_lat_lon_m3, (1, 2, 0))

    collision_freq = np.asarray(base_grid.collision_freq, dtype=float)
    if (
        base_grid.electron_temp_k is not None
        and base_grid.ion_temp_k is not None
        and base_grid.neutral_species_cm3 is not None
    ):
        collision_freq = effective_collision_frequency(
            base_grid.electron_temp_k,
            base_grid.ion_temp_k,
            perturbed_density_lat_lon_alt_m3,
            base_grid.neutral_species_cm3,
        )

    return IonosphereGrid(
        latitudes_deg=base_grid.latitudes_deg,
        longitudes_deg=base_grid.longitudes_deg,
        altitudes_km=base_grid.altitudes_km,
        iono_en_grid=perturbed_density_lat_lon_alt_m3 / 1e6,
        iono_en_grid_5=perturbed_density_lat_lon_alt_m3 / 1e6,
        collision_freq=collision_freq,
        iono_grid_parms=base_grid.iono_grid_parms,
        Bx=base_grid.Bx,
        By=base_grid.By,
        Bz=base_grid.Bz,
        geomag_grid_parms=base_grid.geomag_grid_parms,
        electron_temp_k=base_grid.electron_temp_k,
        ion_temp_k=base_grid.ion_temp_k,
        neutral_temp_k=base_grid.neutral_temp_k,
        neutral_species_cm3=base_grid.neutral_species_cm3,
        metadata=dict(base_grid.metadata),
    )


def _trace_fan_for_case_time(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid: IonosphereGrid,
    fan_elevations_deg: np.ndarray,
    fan_bearings_deg: np.ndarray,
    frequencies_mhz: np.ndarray,
    ox_mode: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nfan = int(fan_elevations_deg.size)
    nfreq = int(frequencies_mhz.size)
    elevs = np.tile(np.asarray(fan_elevations_deg, dtype=float), nfreq)
    bears = np.tile(np.asarray(fan_bearings_deg, dtype=float), nfreq)
    freqs = np.repeat(np.asarray(frequencies_mhz, dtype=float), nfan)

    group_range_km = np.full((nfreq, nfan), np.nan, dtype=float)
    miss_m = np.full((nfreq, nfan), np.nan, dtype=float)
    absorption_db = np.full((nfreq, nfan), np.nan, dtype=float)

    prepared = tracer.prepare_ray_state_vector_batch(
        tx=tx,
        grid=grid,
        elevations_deg=elevs,
        bearings_deg=bears,
        freqs_mhz=freqs,
        ox_mode=ox_mode,
    )
    if prepared is None:
        return group_range_km, miss_m, absorption_db

    tx_local, state_vector, valid_mask = prepared
    if not np.any(valid_mask):
        return group_range_km, miss_m, absorption_db

    rays = tracer.trace_state_vector_batch(
        tx=tx_local,
        elevations_deg=elevs[valid_mask],
        bearings_deg=bears[valid_mask],
        freqs_mhz=freqs[valid_mask],
        grid=grid,
        state_vector=state_vector,
        ox_mode=ox_mode,
        nhops=1,
    )
    for flat_index, ray in zip(np.flatnonzero(valid_mask), rays):
        distance_m, group_path_km, total_absorption_db = _return_leg_distance(ray.path, tx, rx)
        if not math.isfinite(distance_m):
            continue
        freq_index = int(flat_index // nfan)
        fan_index = int(flat_index % nfan)
        group_range_km[freq_index, fan_index] = float(group_path_km) if group_path_km is not None else np.nan
        miss_m[freq_index, fan_index] = float(distance_m)
        absorption_db[freq_index, fan_index] = float(total_absorption_db) if total_absorption_db is not None else 0.0
    return group_range_km, miss_m, absorption_db


def _canonical_launch_angles(elevation_deg: float, bearing_deg: float) -> tuple[float, float]:
    """Map equivalent launch angles to elevation within [-90, 90] degrees."""
    elevation_rad = math.radians(elevation_deg)
    if math.cos(elevation_rad) < 0.0:
        bearing_deg += 180.0
    return math.degrees(math.asin(math.sin(elevation_rad))), bearing_deg % 360.0


def _trace_candidate_topside(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid: IonosphereGrid,
    elevation_deg: float,
    bearing_deg: float,
    frequency_mhz: float,
    ox_mode: int,
    transmitter_state: tuple[GeoPoint, float, float, float, float] | None = None,
) -> tuple[object, float, float | None, float | None, float]:
    # Nelder-Mead can step through the nadir pole. Reflect such angles back
    # into the physical elevation range, rotating the bearing by 180 degrees.
    elevation_deg, bearing_deg = _canonical_launch_angles(elevation_deg, bearing_deg)
    prepared = tracer.prepare_ray_state_vector_batch(
        tx=tx,
        grid=grid,
        elevations_deg=[elevation_deg],
        bearings_deg=[bearing_deg],
        freqs_mhz=[frequency_mhz],
        ox_mode=ox_mode,
        transmitter_state=transmitter_state,
    )
    if prepared is None:
        return None, math.inf, None, None, float("nan")
    tx_local, state_vector, valid_mask = prepared
    if valid_mask.size != 1 or not bool(valid_mask[0]):
        return None, math.inf, None, None, float("nan")
    ray = tracer.trace_state_vector_batch(
        tx=tx_local,
        elevations_deg=[elevation_deg],
        bearings_deg=[bearing_deg],
        freqs_mhz=[frequency_mhz],
        grid=grid,
        state_vector=state_vector,
        ox_mode=ox_mode,
        nhops=1,
    )[0]
    miss_m, group_path_km, absorption_db = _return_leg_distance(ray.path, tx, rx)
    doppler_hz = _doppler_from_summary(ray)
    return ray, miss_m, group_path_km, absorption_db, doppler_hz


def _follow_neighbor_return(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid: IonosphereGrid,
    frequency_mhz: float,
    ox_mode: int,
    start_elevation_deg: float,
    start_bearing_deg: float,
    homing_tolerance_m: float,
) -> tuple[object, float, float | None, float | None, float]:
    cache: dict[tuple[int, int], tuple[object, float, float | None, float | None, float]] = {}
    transmitter_state = tracer.prepare_transmitter_state(tx=tx, grid=grid)

    def evaluate(elevation_deg: float, bearing_deg: float) -> tuple[object, float, float | None, float | None, float]:
        elevation_deg, bearing_deg = _canonical_launch_angles(elevation_deg, bearing_deg)
        key = (int(round(elevation_deg * 1e5)), int(round(bearing_deg * 1e5)))
        if key in cache:
            return cache[key]
        result = _trace_candidate_topside(
            tracer,
            tx=tx,
            rx=rx,
            grid=grid,
            elevation_deg=elevation_deg,
            bearing_deg=bearing_deg,
            frequency_mhz=frequency_mhz,
            ox_mode=ox_mode,
            transmitter_state=transmitter_state,
        )
        cache[key] = result
        return result

    ray, miss_m, group_path_km, absorption_db, doppler_hz = evaluate(start_elevation_deg, start_bearing_deg)
    if math.isfinite(miss_m) and miss_m <= homing_tolerance_m:
        return ray, miss_m, group_path_km, absorption_db, doppler_hz

    def objective(params: np.ndarray) -> float:
        _ray, miss_local_m, _group, _absorption, _doppler = evaluate(float(params[0]), float(params[1]))
        return miss_local_m

    result = minimize(
        objective,
        np.array([start_elevation_deg, start_bearing_deg], dtype=float),
        method="Nelder-Mead",
        options={"fatol": homing_tolerance_m, "xatol": 0.01, "maxfev": 50, "disp": False},
    )
    return evaluate(float(result.x[0]), float(result.x[1]))


def _home_frequency_returns(
    tracer: PointToPointRayTracer,
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid: IonosphereGrid,
    fan_elevations_deg: np.ndarray,
    fan_bearings_deg: np.ndarray,
    frequency_mhz: float,
    ox_mode: int,
    config: TopsideInverseConfig,
    optimizer_method: str = "Nelder-Mead",
) -> tuple[HomedRayReturn, ...]:
    if optimizer_method not in ("Nelder-Mead", "Powell"):
        raise ValueError("optimizer_method must be Nelder-Mead or Powell")
    elevs = np.asarray(fan_elevations_deg, dtype=float)
    bears = np.asarray(fan_bearings_deg, dtype=float)
    transmitter_state = tracer.prepare_transmitter_state(tx=tx, grid=grid)
    prepared = tracer.prepare_ray_state_vector_batch(
        tx=tx,
        grid=grid,
        elevations_deg=elevs,
        bearings_deg=bears,
        freqs_mhz=np.full(elevs.shape, float(frequency_mhz), dtype=float),
        ox_mode=ox_mode,
        transmitter_state=transmitter_state,
    )
    if prepared is None:
        return ()
    tx_local, state_vector, valid_mask = prepared
    rays = tracer.trace_state_vector_batch(
        tx=tx_local,
        elevations_deg=elevs[valid_mask],
        bearings_deg=bears[valid_mask],
        freqs_mhz=np.full(int(np.count_nonzero(valid_mask)), float(frequency_mhz), dtype=float),
        grid=grid,
        state_vector=state_vector,
        ox_mode=ox_mode,
        nhops=1,
    )
    miss_flat = np.full(elevs.shape, np.inf, dtype=float)
    for flat_index, ray in zip(np.flatnonzero(valid_mask), rays):
        miss_m, _group_path_km, _absorption_db = _return_leg_distance(ray.path, tx, rx)
        miss_flat[int(flat_index)] = miss_m

    unique_elevs = np.unique(elevs)
    unique_bears = np.unique(bears)
    if unique_elevs.size * unique_bears.size != elevs.size:
        return ()
    miss_grid = miss_flat.reshape(unique_bears.size, unique_elevs.size)
    local_minima = _local_minimum_mask(miss_grid)
    seed_indices = np.flatnonzero(local_minima.ravel())
    if seed_indices.size == 0:
        finite = np.flatnonzero(np.isfinite(miss_flat))
        seed_indices = finite[: config.seed_max_candidates_per_frequency]
    else:
        seed_indices = seed_indices[np.argsort(miss_flat[seed_indices])]
        seed_indices = seed_indices[: config.seed_max_candidates_per_frequency]

    cache: dict[tuple[int, int], tuple[object, float, float | None, float | None, float]] = {}

    def evaluate(elevation_deg: float, bearing_deg: float) -> tuple[object, float, float | None, float | None, float]:
        elevation_deg, bearing_deg = _canonical_launch_angles(elevation_deg, bearing_deg)
        key = (int(round(elevation_deg * 1e5)), int(round(bearing_deg * 1e5)))
        cached = cache.get(key)
        if cached is not None:
            return cached
        result = _trace_candidate_topside(
            tracer,
            tx=tx,
            rx=rx,
            grid=grid,
            elevation_deg=elevation_deg,
            bearing_deg=bearing_deg,
            frequency_mhz=frequency_mhz,
            ox_mode=ox_mode,
            transmitter_state=transmitter_state,
        )
        cache[key] = result
        return result

    seed_points: list[tuple[float, float]] = [
        (float(elevs[int(flat_index)]), float(bears[int(flat_index)])) for flat_index in seed_indices
    ]
    if tx.alt_km > 400.0 and rx.alt_km > 400.0:
        los_bearing_deg, los_elevation_deg, _ = _line_of_sight_angles(tx, rx)
        seed_points.insert(0, (float(los_elevation_deg), float(los_bearing_deg)))

    vertical_fan = bool(np.all(elevs < 0.0))
    near_nadir_limit = float(unique_elevs[min(1, unique_elevs.size - 1)])

    def home_seed(start_elev: float, start_bear: float, seed_method: str) -> HomedRayReturn | None:
        use_nadir_coordinates = seed_method == "Powell" and vertical_fan

        def angles_from_params(params: np.ndarray) -> tuple[float, float]:
            if not use_nadir_coordinates:
                return float(params[0]), float(params[1])
            nadir_offset = math.hypot(float(params[0]), float(params[1]))
            return -90.0 + nadir_offset, math.degrees(math.atan2(float(params[0]), float(params[1]))) % 360.0

        def objective(params: np.ndarray) -> float:
            if use_nadir_coordinates and math.hypot(float(params[0]), float(params[1])) > 90.0:
                return 1e9
            elevation, bearing = angles_from_params(params)
            _ray, miss_m, _group, _absorption, _doppler = evaluate(elevation, bearing)
            return miss_m

        if seed_method == "Powell":
            options = {"ftol": 0.01, "xtol": 0.01, "maxfev": 80, "disp": False}
            nadir_offset = 90.0 + start_elev
            bearing_rad = math.radians(start_bear)
            start_params = (np.array([nadir_offset * math.sin(bearing_rad),
                                      nadir_offset * math.cos(bearing_rad)]) if use_nadir_coordinates
                            else np.array([start_elev, start_bear], dtype=float))
        else:
            options = {"fatol": config.homing_tolerance_m, "xatol": 0.01, "maxfev": 80, "disp": False}
            start_params = np.array([start_elev, start_bear], dtype=float)
        result = minimize(
            objective,
            start_params,
            method=seed_method,
            options=options,
        )
        result_elev, result_bear = angles_from_params(result.x)
        ray, miss_m, group_path_km, absorption_db, doppler_hz = evaluate(result_elev, result_bear)
        if (use_nadir_coordinates and 10.0 < miss_m < 5_000.0
                and group_path_km is not None and group_path_km >= config.range_min_km):
            def polish_objective(params: np.ndarray) -> float:
                return evaluate(float(params[0]), float(params[1]))[1]

            polished = minimize(
                polish_objective, np.array([result_elev, result_bear]),
                method="Nelder-Mead",
                options={"fatol": 10.0, "xatol": 0.01, "maxfev": 40, "disp": False},
            )
            polished_ray = evaluate(float(polished.x[0]), float(polished.x[1]))
            if polished_ray[1] < miss_m:
                ray, miss_m, group_path_km, absorption_db, doppler_hz = polished_ray
        if ray is None or not math.isfinite(miss_m) or group_path_km is None or not math.isfinite(group_path_km):
            return None
        if miss_m > config.homing_tolerance_m:
            return None
        return HomedRayReturn(
            frequency_mhz=float(frequency_mhz),
            ox_mode=int(ox_mode),
            ray=ray,
            miss_m=float(miss_m),
            group_range_km=float(group_path_km),
            absorption_db=0.0 if absorption_db is None or not math.isfinite(absorption_db) else float(absorption_db),
            doppler_hz=float(doppler_hz) if math.isfinite(doppler_hz) else 0.0,
        )

    homed: list[HomedRayReturn] = []
    for start_elev, start_bear in seed_points:
        seed_method = ("Nelder-Mead" if optimizer_method == "Powell" and vertical_fan
                       and start_elev <= near_nadir_limit + 1e-9 else optimizer_method)
        solution = home_seed(start_elev, start_bear, seed_method)
        if solution is not None:
            homed.append(solution)

    return _deduplicate_homed_returns(homed, config.homed_max_returns_per_frequency)


def _deduplicate_homed_returns(
    returns: Sequence[HomedRayReturn], max_returns: int,
) -> tuple[HomedRayReturn, ...]:
    """Keep distinct physical launch directions and virtual ranges."""
    chosen: list[HomedRayReturn] = []
    chosen_directions: list[np.ndarray] = []
    angular_cosine = math.cos(math.radians(0.02))
    for solution in sorted(returns, key=lambda item: (item.group_range_km, item.miss_m)):
        elevation = math.radians(float(solution.ray.path["initial_elev"]))
        bearing = math.radians(float(solution.ray.path["initial_bearing"]))
        direction = np.array([
            math.cos(elevation) * math.sin(bearing),
            math.cos(elevation) * math.cos(bearing),
            math.sin(elevation),
        ])
        if any(
            abs(solution.group_range_km - prior.group_range_km) < 1.0
            and float(np.dot(direction, prior_direction)) >= angular_cosine
            for prior, prior_direction in zip(chosen, chosen_directions)
        ):
            continue
        chosen.append(solution)
        chosen_directions.append(direction)
        if len(chosen) >= max_returns:
            break
    return tuple(chosen)


def _continue_homed_returns(
    tracer: PointToPointRayTracer,
    previous: Sequence[HomedRayReturn],
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid: IonosphereGrid,
    frequency_mhz: float,
    ox_mode: int,
    config: TopsideInverseConfig,
    range_min_km: float,
) -> tuple[HomedRayReturn, ...]:
    """Home the previous frequency's paths at an adjacent frequency."""
    followed: list[HomedRayReturn] = []
    for solution in previous:
        ray, miss_m, group_range_km, absorption_db, doppler_hz = _follow_neighbor_return(
            tracer,
            tx=tx,
            rx=rx,
            grid=grid,
            frequency_mhz=frequency_mhz,
            ox_mode=ox_mode,
            start_elevation_deg=float(solution.ray.path["initial_elev"]),
            start_bearing_deg=float(solution.ray.path["initial_bearing"]),
            homing_tolerance_m=config.homing_tolerance_m,
        )
        if (ray is None or not math.isfinite(miss_m) or miss_m > config.homing_tolerance_m
                or group_range_km is None or not math.isfinite(group_range_km)
                or group_range_km < range_min_km):
            continue
        followed.append(HomedRayReturn(
            frequency_mhz=frequency_mhz,
            ox_mode=ox_mode,
            ray=ray,
            miss_m=float(miss_m),
            group_range_km=float(group_range_km),
            absorption_db=0.0 if absorption_db is None or not math.isfinite(absorption_db) else float(absorption_db),
            doppler_hz=float(doppler_hz) if math.isfinite(doppler_hz) else 0.0,
        ))
    return _deduplicate_homed_returns(followed, config.homed_max_returns_per_frequency)


def _cached_topside_tracer() -> PointToPointRayTracer:
    return PointToPointRayTracer(cache_native_grid=True)


def home_frequency_sweep_adaptive(
    *,
    tx: GeoPoint,
    rx: GeoPoint,
    grid: IonosphereGrid,
    fan_elevations_deg: np.ndarray,
    fan_bearings_deg: np.ndarray,
    frequencies_mhz: Sequence[float],
    ox_mode: int,
    config: TopsideInverseConfig,
    range_min_km: float = 150.0,
    anchor_stride: int = 5,
    anchor_block_size: int = 10,
    anchor_elevation_stride: int = 2,
    anchor_optimizer: str = "Powell",
    tracer_factory: Callable[[], PointToPointRayTracer] = _cached_topside_tracer,
) -> tuple[tuple[HomedRayReturn, ...], ...]:
    """Trace a frequency sweep with periodic fan anchors and neighboring-ray homing.

    Every ``anchor_stride`` steps, each ``anchor_block_size`` block end, and
    the final step search an elevation-subset fan. The first two rows nearest
    nadir are always retained. Powell anchors retain Nelder-Mead for those
    rows. Between anchors, accepted returns are followed forward and backward.
    A failed continuation triggers an anchor-fan search.
    Only returns at or above ``range_min_km`` are retained.
    """
    frequencies = np.asarray(frequencies_mhz, dtype=float)
    if frequencies.ndim != 1 or np.any(np.diff(frequencies) <= 0.0):
        raise ValueError("frequencies_mhz must be a strictly increasing vector")
    if anchor_stride < 1:
        raise ValueError("anchor_stride must be positive")
    if anchor_block_size < 1:
        raise ValueError("anchor_block_size must be positive")
    if anchor_elevation_stride < 1:
        raise ValueError("anchor_elevation_stride must be positive")
    if frequencies.size == 0:
        return ()

    anchor_elevations = np.asarray(fan_elevations_deg, dtype=float)
    anchor_bearings = np.asarray(fan_bearings_deg, dtype=float)
    if anchor_elevation_stride > 1:
        unique_elevations = np.unique(anchor_elevations)
        if unique_elevations.size * np.unique(anchor_bearings).size == anchor_elevations.size:
            retained = np.zeros(unique_elevations.size, dtype=bool)
            retained[:2] = True
            retained[2::anchor_elevation_stride] = True
            retained[-1] = True
            mask = np.isin(anchor_elevations, unique_elevations[retained])
            anchor_elevations = anchor_elevations[mask]
            anchor_bearings = anchor_bearings[mask]

    results: list[tuple[HomedRayReturn, ...]] = []
    previous_anchor = 0

    def global_search(frequency_mhz: float) -> tuple[HomedRayReturn, ...]:
        found = _home_frequency_returns(
            tracer_factory(), tx=tx, rx=rx, grid=grid,
            fan_elevations_deg=anchor_elevations,
            fan_bearings_deg=anchor_bearings,
            frequency_mhz=frequency_mhz, ox_mode=ox_mode, config=config,
            optimizer_method=anchor_optimizer,
        )
        return tuple(solution for solution in found if solution.group_range_km >= range_min_km)

    for frequency_index, frequency_mhz in enumerate(frequencies):
        frequency_mhz = float(frequency_mhz)
        is_anchor = (frequency_index % anchor_stride == 0
                     or frequency_index % anchor_block_size == anchor_block_size - 1
                     or frequency_index == frequencies.size - 1)
        if is_anchor:
            current = global_search(frequency_mhz)
        else:
            current = _continue_homed_returns(
                tracer_factory(), results[-1], tx=tx, rx=rx, grid=grid,
                frequency_mhz=frequency_mhz, ox_mode=ox_mode,
                config=config, range_min_km=range_min_km,
            )
            if not current and results[-1]:
                current = global_search(frequency_mhz)
        results.append(current)

        if is_anchor and frequency_index > previous_anchor:
            for back_index in range(frequency_index - 1, previous_anchor, -1):
                backward = _continue_homed_returns(
                    tracer_factory(), results[back_index + 1], tx=tx, rx=rx,
                    grid=grid, frequency_mhz=float(frequencies[back_index]),
                    ox_mode=ox_mode, config=config, range_min_km=range_min_km,
                )
                results[back_index] = _deduplicate_homed_returns(
                    (*results[back_index], *backward),
                    config.homed_max_returns_per_frequency,
                )
                if not results[back_index] and results[back_index + 1]:
                    results[back_index] = global_search(float(frequencies[back_index]))
            previous_anchor = frequency_index
    return tuple(results)


def _estimate_return_doppler_hz(
    tracer: PointToPointRayTracer,
    *,
    case: TopsideCase,
    center_return: HomedRayReturn,
    grid: IonosphereGrid,
    config: TopsideInverseConfig,
) -> float:
    dt_seconds = 0.5 * (case.times_utc[2] - case.times_utc[0]).total_seconds()
    if dt_seconds <= 0.0:
        return 0.0
    launch_elev = float(center_return.ray.path.get("initial_elev", np.nan))
    launch_bear = float(center_return.ray.path.get("initial_bearing", np.nan))
    if not (math.isfinite(launch_elev) and math.isfinite(launch_bear)):
        return 0.0

    minus = _follow_neighbor_return(
        tracer,
        tx=case.tx_points[0],
        rx=case.rx_points[0],
        grid=grid,
        frequency_mhz=center_return.frequency_mhz,
        ox_mode=center_return.ox_mode,
        start_elevation_deg=launch_elev,
        start_bearing_deg=launch_bear,
        homing_tolerance_m=config.homing_tolerance_m,
    )
    plus = _follow_neighbor_return(
        tracer,
        tx=case.tx_points[2],
        rx=case.rx_points[2],
        grid=grid,
        frequency_mhz=center_return.frequency_mhz,
        ox_mode=center_return.ox_mode,
        start_elevation_deg=launch_elev,
        start_bearing_deg=launch_bear,
        homing_tolerance_m=config.homing_tolerance_m,
    )
    group_minus_km = minus[2]
    group_plus_km = plus[2]
    if group_minus_km is None or group_plus_km is None:
        return 0.0
    if not (math.isfinite(group_minus_km) and math.isfinite(group_plus_km)):
        return 0.0
    range_rate_mps = (float(group_plus_km) - float(group_minus_km)) * 1000.0 / (2.0 * dt_seconds)
    return -float(center_return.frequency_mhz) * 1e6 * range_rate_mps / C_M_PER_S


def _accumulate_image(
    image: np.ndarray,
    *,
    returns: Sequence[HomedRayReturn],
    frequencies_mhz: np.ndarray,
    range_edges_km: np.ndarray,
) -> None:
    freq_lookup = {float(freq): index for index, freq in enumerate(np.asarray(frequencies_mhz, dtype=float))}
    for solution in returns:
        freq_index = freq_lookup.get(float(solution.frequency_mhz))
        if freq_index is None:
            continue
        range_bin = int(np.digitize(solution.group_range_km, range_edges_km) - 1)
        if range_bin < 0 or range_bin >= image.shape[1]:
            continue
        weight = 10.0 ** (-max(solution.absorption_db, 0.0) / 10.0)
        image[freq_index, range_bin] += weight


def _accumulate_doppler(
    numerator: np.ndarray,
    denominator: np.ndarray,
    *,
    returns: Sequence[HomedRayReturn],
    frequencies_mhz: np.ndarray,
    range_edges_km: np.ndarray,
) -> None:
    freq_lookup = {float(freq): index for index, freq in enumerate(np.asarray(frequencies_mhz, dtype=float))}
    for solution in returns:
        freq_index = freq_lookup.get(float(solution.frequency_mhz))
        if freq_index is None:
            continue
        range_bin = int(np.digitize(solution.group_range_km, range_edges_km) - 1)
        if range_bin < 0 or range_bin >= numerator.shape[1]:
            continue
        weight = 10.0 ** (-max(solution.absorption_db, 0.0) / 10.0)
        numerator[freq_index, range_bin] += weight * solution.doppler_hz
        denominator[freq_index, range_bin] += weight


def _render_case_observables(
    problem: InverseProblem,
    case: TopsideCase,
    grid: IonosphereGrid,
    tracer: PointToPointRayTracer | None,
) -> CaseObservables:
    total_shape = (problem.frequencies_mhz.size, problem.range_centers_km.size)
    o_image = np.zeros(total_shape, dtype=float)
    x_image = np.zeros(total_shape, dtype=float)
    doppler_num = np.zeros(total_shape, dtype=float)
    doppler_den = np.zeros(total_shape, dtype=float)
    center_tx = case.tx_points[1]
    center_rx = case.rx_points[1]
    for ox_mode in (1, -1):
        for frequency_mhz in problem.frequencies_mhz:
            local_tracer = tracer if tracer is not None else PointToPointRayTracer()
            frequency_returns = _home_frequency_returns(
                local_tracer,
                tx=center_tx,
                rx=center_rx,
                grid=grid,
                fan_elevations_deg=case.fan_elevations_deg,
                fan_bearings_deg=case.fan_bearings_deg,
                frequency_mhz=float(frequency_mhz),
                ox_mode=ox_mode,
                config=problem.config,
            )
            if ox_mode == 1:
                _accumulate_image(
                    o_image,
                    returns=frequency_returns,
                    frequencies_mhz=problem.frequencies_mhz,
                    range_edges_km=problem.range_edges_km,
                )
            else:
                _accumulate_image(
                    x_image,
                    returns=frequency_returns,
                    frequencies_mhz=problem.frequencies_mhz,
                    range_edges_km=problem.range_edges_km,
                )
            frequency_returns_with_doppler = [
                replace(
                    solution,
                    doppler_hz=_estimate_return_doppler_hz(
                        local_tracer,
                        case=case,
                        center_return=solution,
                        grid=grid,
                        config=problem.config,
                    ),
                )
                for solution in frequency_returns
            ]
            _accumulate_doppler(
                doppler_num,
                doppler_den,
                returns=frequency_returns_with_doppler,
                frequencies_mhz=problem.frequencies_mhz,
                range_edges_km=problem.range_edges_km,
            )
            del frequency_returns_with_doppler
            del frequency_returns
            if tracer is None:
                del local_tracer
            gc.collect()

    sigma = (problem.config.blur_sigma_frequency_bins, problem.config.blur_sigma_range_bins)
    o_image = gaussian_filter(o_image, sigma=sigma, mode="nearest")
    x_image = gaussian_filter(x_image, sigma=sigma, mode="nearest")
    doppler_num = gaussian_filter(doppler_num, sigma=sigma, mode="nearest")
    doppler_den = gaussian_filter(doppler_den, sigma=sigma, mode="nearest")

    total_image = o_image + x_image
    total_scale = max(float(np.max(total_image)), 1e-9)
    o_norm = o_image / total_scale
    x_norm = x_image / total_scale
    total_norm = total_image / total_scale
    ox_split = (o_image - x_image) / np.maximum(total_image, 1e-6)
    doppler_image_hz = np.divide(doppler_num, doppler_den, out=np.zeros_like(doppler_num), where=doppler_den > 1e-9)

    return CaseObservables(
        name=case.name,
        kind=case.kind,
        total_image=total_norm,
        o_image=o_norm,
        x_image=x_norm,
        ox_split_image=ox_split,
        doppler_image_hz=doppler_image_hz,
        doppler_weight=np.clip(total_norm, 0.0, 1.0),
    )


def simulate_dataset(problem: InverseProblem, params: IonosphereFitParams) -> SyntheticDataset:
    cases = tuple(
        _render_case_observables(
            problem,
            case,
            _apply_fit_params_to_grid(problem, background_grid, params),
            None,
        )
        for case, background_grid in zip(problem.cases, problem.background_grids)
    )
    return SyntheticDataset(
        frequencies_mhz=problem.frequencies_mhz,
        range_edges_km=problem.range_edges_km,
        range_centers_km=problem.range_centers_km,
        cases=cases,
    )


def dataset_cost(observed: SyntheticDataset, predicted: SyntheticDataset) -> float:
    case_costs: list[float] = []
    doppler_scale_hz = 50.0
    for observed_case, predicted_case in zip(observed.cases, predicted.cases):
        image_cost = float(np.mean((predicted_case.total_image - observed_case.total_image) ** 2))
        split_weight = np.maximum(observed_case.total_image, 0.0)
        if float(np.sum(split_weight)) > 0.0:
            split_cost = float(
                np.sum(((predicted_case.ox_split_image - observed_case.ox_split_image) ** 2) * split_weight)
                / np.sum(split_weight)
            )
        else:
            split_cost = 0.0
        doppler_weight = np.maximum(observed_case.doppler_weight, 0.0)
        if float(np.sum(doppler_weight)) > 0.0:
            doppler_cost = float(
                np.sum(
                    (((predicted_case.doppler_image_hz - observed_case.doppler_image_hz) / doppler_scale_hz) ** 2)
                    * doppler_weight
                )
                / np.sum(doppler_weight)
            )
        else:
            doppler_cost = 0.0
        case_costs.append(image_cost + 0.7 * doppler_cost + 0.5 * split_cost)
    return float(np.mean(case_costs)) if case_costs else float("inf")


def solve_inverse_problem(problem: InverseProblem, observed: SyntheticDataset) -> FitResult:
    truth_prediction = simulate_dataset(problem, problem.config.truth_params)
    truth_cost = dataset_cost(observed, truth_prediction)
    del truth_prediction
    gc.collect()
    return solve_inverse_problem_from_observed(problem, observed, truth_cost=truth_cost)


def solve_inverse_problem_from_observed(
    problem: InverseProblem,
    observed: SyntheticDataset,
    *,
    truth_cost: float,
) -> FitResult:
    parameter_names = tuple(problem.config.fit_parameter_names)
    if problem.config.solver_maxiter <= 0:
        fitted_params = problem.config.initial_params
        fitted_prediction = simulate_dataset(problem, fitted_params)
        fitted_cost = dataset_cost(observed, fitted_prediction)
        del fitted_prediction
        gc.collect()
        return FitResult(
            truth_params=problem.config.truth_params,
            fitted_params=fitted_params,
            truth_cost=float(truth_cost),
            fitted_cost=float(fitted_cost),
            solver_success=True,
            solver_message="solver skipped because solver_maxiter <= 0; using initial_params",
            iterations=0,
            evaluations=1,
        )
    parameter_bounds = _fit_parameter_bounds()
    bounds = [parameter_bounds[name] for name in parameter_names]

    def objective(vector: np.ndarray) -> float:
        params = IonosphereFitParams.from_vector(vector, parameter_names, base=problem.config.initial_params)
        predicted = simulate_dataset(problem, params)
        cost = dataset_cost(observed, predicted)
        del predicted
        gc.collect()
        return cost

    result = differential_evolution(
        objective,
        bounds,
        maxiter=max(problem.config.solver_maxiter, 0),
        popsize=max(problem.config.solver_popsize, 4),
        seed=int(problem.config.solver_seed),
        polish=False,
        updating="deferred",
        workers=1,
        disp=False,
    )
    fitted_params = IonosphereFitParams.from_vector(result.x, parameter_names, base=problem.config.initial_params)
    fitted_cost = float(result.fun)
    return FitResult(
        truth_params=problem.config.truth_params,
        fitted_params=fitted_params,
        truth_cost=float(truth_cost),
        fitted_cost=fitted_cost,
        solver_success=bool(result.success),
        solver_message=str(result.message),
        iterations=int(result.nit),
        evaluations=int(result.nfev),
    )


def run_inverse_fit_only(
    config: TopsideInverseConfig | None = None,
    *,
    fit_case_names: Sequence[str] | None = None,
    keep_global_background: bool = False,
) -> FitResult:
    config = config or TopsideInverseConfig()
    problem = _subset_problem_cases(build_inverse_problem(config), fit_case_names)
    if not keep_global_background and problem.global_background_grid is not None:
        problem = replace(problem, global_background_grid=None)
        gc.collect()
    observed = simulate_dataset(problem, config.truth_params)
    return solve_inverse_problem_from_observed(problem, observed, truth_cost=0.0)


def run_inverse_demo(
    config: TopsideInverseConfig | None = None,
    *,
    fit_case_names: Sequence[str] | None = None,
    keep_global_background: bool = False,
) -> InverseDemoResult:
    config = config or TopsideInverseConfig()
    problem = _subset_problem_cases(build_inverse_problem(config), fit_case_names)
    if not keep_global_background and problem.global_background_grid is not None:
        problem = replace(problem, global_background_grid=None)
        gc.collect()
    observed = simulate_dataset(problem, config.truth_params)
    fit = solve_inverse_problem_from_observed(problem, observed, truth_cost=0.0)
    fitted = simulate_dataset(problem, fit.fitted_params)
    return InverseDemoResult(
        problem=problem,
        observed=observed,
        fitted=fitted,
        fit=fit,
    )


def _selected_case_index(demo: InverseDemoResult, case_name: str | None = None) -> int:
    if case_name is None:
        return max(
            range(len(demo.observed.cases)),
            key=lambda index: _plot_case_score(demo.observed.cases[index]),
        )
    matches = [index for index, case in enumerate(demo.problem.cases) if case.name == case_name]
    if not matches:
        raise ValueError(f"unknown case name {case_name!r}")
    return matches[0]


def _plot_case_score(case: CaseObservables) -> float:
    support = np.asarray(case.doppler_weight, dtype=float)
    doppler = np.asarray(case.doppler_image_hz, dtype=float)
    finite_support = support[np.isfinite(support)]
    if finite_support.size == 0:
        return float(np.nanmax(np.asarray(case.total_image, dtype=float)))
    support_threshold = max(1e-6, 0.05 * float(np.nanmax(finite_support)))
    visible = np.isfinite(doppler) & np.isfinite(support) & (support >= support_threshold)
    if np.any(visible):
        return float(np.nanmax(np.abs(doppler[visible])) * np.sum(support[visible]))
    return float(np.nanmax(np.asarray(case.total_image, dtype=float)))


def _masked_visible_field(
    field: np.ndarray,
    support: np.ndarray,
    *,
    threshold_fraction: float = 0.05,
) -> tuple[np.ma.MaskedArray, np.ndarray]:
    field = np.asarray(field, dtype=float)
    support = np.asarray(support, dtype=float)
    finite_support = support[np.isfinite(support)]
    if finite_support.size == 0:
        visible = np.zeros(field.shape, dtype=bool)
    else:
        support_threshold = max(1e-6, threshold_fraction * float(np.nanmax(finite_support)))
        visible = np.isfinite(field) & np.isfinite(support) & (support >= support_threshold)
    return np.ma.array(field, mask=~visible), field[visible]


def _max_observed_range_km(
    image: np.ndarray,
    range_centers_km: np.ndarray,
    range_edges_km: np.ndarray,
    *,
    threshold_fraction: float = 0.05,
) -> float:
    field = np.asarray(image, dtype=float)
    finite = field[np.isfinite(field)]
    if finite.size == 0:
        return float(range_edges_km[-1])
    threshold = max(1e-6, threshold_fraction * float(np.max(finite)))
    col_support = np.any(np.isfinite(field) & (field >= threshold), axis=0)
    if not np.any(col_support):
        return float(range_edges_km[-1])
    last_index = int(np.flatnonzero(col_support)[-1])
    last_index = min(last_index + 1, len(range_edges_km) - 1)
    return float(range_edges_km[last_index])


def _observed_support_extents(
    image: np.ndarray,
    frequencies_mhz: np.ndarray,
    range_edges_km: np.ndarray,
    *,
    threshold_fraction: float = 0.05,
) -> tuple[float, float, float, float]:
    field = np.asarray(image, dtype=float)
    finite = field[np.isfinite(field)]
    if finite.size == 0:
        return (
            float(frequencies_mhz[0]),
            float(frequencies_mhz[-1]),
            float(range_edges_km[0]),
            float(range_edges_km[-1]),
        )
    threshold = max(1e-6, threshold_fraction * float(np.max(finite)))
    support = np.isfinite(field) & (field >= threshold)
    if not np.any(support):
        return (
            float(frequencies_mhz[0]),
            float(frequencies_mhz[-1]),
            float(range_edges_km[0]),
            float(range_edges_km[-1]),
        )
    freq_rows = np.flatnonzero(np.any(support, axis=1))
    range_cols = np.flatnonzero(np.any(support, axis=0))
    fmin = float(frequencies_mhz[int(freq_rows[0])])
    fmax = float(frequencies_mhz[int(freq_rows[-1])])
    rmin = float(range_edges_km[int(range_cols[0])])
    rmax = float(range_edges_km[min(int(range_cols[-1]) + 1, len(range_edges_km) - 1)])
    return fmin, fmax, rmin, rmax


def _dense_plot_xlim_mhz(
    frequencies_mhz: np.ndarray,
    range_centers_km: np.ndarray,
    *power_fields: np.ndarray,
    threshold_fraction: float = 0.05,
) -> tuple[float, float]:
    frequencies = np.asarray(frequencies_mhz, dtype=float)
    if frequencies.size == 0:
        return 0.0, 1.0
    if frequencies.size == 1:
        return float(frequencies[0] - 0.25), float(frequencies[0] + 0.25)
    step_mhz = float(np.median(np.diff(frequencies)))
    xmin = float(frequencies[0] - 0.5 * step_mhz)
    xmax_full = float(frequencies[-1] + 0.5 * step_mhz)

    combined = np.maximum.reduce([np.asarray(field, dtype=float) for field in power_fields])
    finite = combined[np.isfinite(combined)]
    if finite.size == 0:
        return xmin, xmax_full
    threshold = max(1e-6, threshold_fraction * float(np.max(finite)))
    support = np.isfinite(combined) & (combined >= threshold)
    if not np.any(support):
        return xmin, xmax_full

    ranges_km = np.asarray(range_centers_km, dtype=float)
    direct_range_km = float(np.min(ranges_km[np.any(support, axis=0)]))
    range_step_km = float(np.median(np.diff(ranges_km))) if ranges_km.size > 1 else 25.0
    refracted_support = support & (ranges_km[None, :] >= direct_range_km + max(75.0, 2.0 * range_step_km))
    if not np.any(refracted_support):
        return xmin, xmax_full

    refracted_freqs = frequencies[np.any(refracted_support, axis=1)]
    xmax = min(xmax_full, float(np.max(refracted_freqs) + 2.5 * step_mhz))
    return xmin, max(xmax, xmin + step_mhz)


def _visible_color_limits(
    field: np.ndarray,
    *,
    x_centers: np.ndarray,
    y_centers: np.ndarray,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
    floor: float | None = None,
) -> tuple[float, float]:
    xmin, xmax = sorted((float(x_limits[0]), float(x_limits[1])))
    ymin, ymax = sorted((float(y_limits[0]), float(y_limits[1])))
    xmask = (np.asarray(x_centers, dtype=float) >= xmin) & (np.asarray(x_centers, dtype=float) <= xmax)
    ymask = (np.asarray(y_centers, dtype=float) >= ymin) & (np.asarray(y_centers, dtype=float) <= ymax)
    if np.any(xmask) and np.any(ymask):
        finite = np.asarray(field, dtype=float)[np.ix_(ymask, xmask)]
    else:
        finite = np.asarray(field, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        vmin = 0.0 if floor is None else float(floor)
        return vmin, vmin + 1.0
    vmin = float(np.min(finite))
    vmax = float(np.max(finite))
    if floor is not None:
        vmin = max(float(floor), vmin)
    if not vmax > vmin:
        vmax = vmin + 1e-6
    return vmin, vmax


def _resample_vertical_display(
    field_alt_x: np.ndarray,
    altitudes_km: np.ndarray,
    *,
    step_km: float,
) -> tuple[np.ndarray, np.ndarray]:
    altitudes = np.asarray(altitudes_km, dtype=float)
    field = np.asarray(field_alt_x, dtype=float)
    if altitudes.size < 2:
        return altitudes, field
    step = max(float(step_km), 1e-3)
    display_alts = np.arange(float(altitudes[0]), float(altitudes[-1]) + 0.5 * step, step, dtype=float)
    resampled = np.empty((display_alts.size, field.shape[1]), dtype=float)
    for ix in range(field.shape[1]):
        resampled[:, ix] = np.interp(
            display_alts,
            altitudes,
            field[:, ix],
            left=float(field[0, ix]),
            right=float(field[-1, ix]),
        )
    return display_alts, resampled


def _case_problem_with_frequencies(
    problem: InverseProblem,
    frequencies_mhz: np.ndarray,
    *,
    range_bin_km: float | None = None,
) -> InverseProblem:
    if range_bin_km is None:
        return replace(problem, frequencies_mhz=np.asarray(frequencies_mhz, dtype=float))
    bin_km = max(float(range_bin_km), 1e-3)
    range_edges_km = np.arange(
        problem.config.range_min_km,
        problem.config.range_max_km + 0.5 * bin_km,
        bin_km,
        dtype=float,
    )
    range_centers_km = 0.5 * (range_edges_km[:-1] + range_edges_km[1:])
    return replace(
        problem,
        frequencies_mhz=np.asarray(frequencies_mhz, dtype=float),
        range_edges_km=range_edges_km,
        range_centers_km=range_centers_km,
    )


def _trace_case_returns(
    problem: InverseProblem,
    case: TopsideCase,
    background_grid: IonosphereGrid,
    params: IonosphereFitParams,
) -> tuple[IonosphereGrid, tuple[HomedRayReturn, ...], tuple[HomedRayReturn, ...]]:
    tracer = PointToPointRayTracer()
    grid = _apply_fit_params_to_grid(problem, background_grid, params)
    center_tx = case.tx_points[1]
    center_rx = case.rx_points[1]
    o_returns: list[HomedRayReturn] = []
    x_returns: list[HomedRayReturn] = []
    for ox_mode in (1, -1):
        mode_returns: list[HomedRayReturn] = []
        for frequency_mhz in problem.frequencies_mhz:
            mode_returns.extend(
                _home_frequency_returns(
                    tracer,
                    tx=center_tx,
                    rx=center_rx,
                    grid=grid,
                    fan_elevations_deg=case.fan_elevations_deg,
                    fan_bearings_deg=case.fan_bearings_deg,
                    frequency_mhz=float(frequency_mhz),
                    ox_mode=ox_mode,
                    config=problem.config,
                )
            )
        if ox_mode == 1:
            o_returns = mode_returns
        else:
            x_returns = mode_returns
    return grid, tuple(o_returns), tuple(x_returns)


def _trace_case_profile_curves(
    problem: InverseProblem,
    case: TopsideCase,
    background_grid: IonosphereGrid,
    params: IonosphereFitParams,
) -> tuple[IonosphereGrid, tuple[RayProfileCurve, ...], tuple[RayProfileCurve, ...]]:
    grid = _apply_fit_params_to_grid(problem, background_grid, params)
    center_tx = case.tx_points[1]
    center_rx = case.rx_points[1]
    o_curves: list[RayProfileCurve] = []
    x_curves: list[RayProfileCurve] = []
    for ox_mode in (1, -1):
        mode_curves = o_curves if ox_mode == 1 else x_curves
        for frequency_mhz in problem.frequencies_mhz:
            tracer = PointToPointRayTracer()
            returns = _home_frequency_returns(
                tracer,
                tx=center_tx,
                rx=center_rx,
                grid=grid,
                fan_elevations_deg=case.fan_elevations_deg,
                fan_bearings_deg=case.fan_bearings_deg,
                frequency_mhz=float(frequency_mhz),
                ox_mode=ox_mode,
                config=problem.config,
            )
            for solution in returns:
                along_km, alt_km = _ray_profile_along_track_km(case, solution.ray)
                if along_km.size == 0 or alt_km.size == 0:
                    continue
                mode_curves.append(
                    RayProfileCurve(
                        frequency_mhz=float(solution.frequency_mhz),
                        ox_mode=int(solution.ox_mode),
                        miss_m=float(solution.miss_m),
                        group_range_km=float(solution.group_range_km),
                        along_km=along_km,
                        alt_km=alt_km,
                    )
                )
            del returns
            del tracer
            gc.collect()
    return grid, tuple(o_curves), tuple(x_curves)


def _render_case_for_params(
    problem: InverseProblem,
    case: TopsideCase,
    background_grid: IonosphereGrid,
    params: IonosphereFitParams,
) -> CaseObservables:
    return _render_case_observables(
        problem,
        case,
        _apply_fit_params_to_grid(problem, background_grid, params),
        None,
    )


def _filter_ionospheric_returns(
    returns: Sequence[HomedRayReturn],
    range_centers_km: np.ndarray,
) -> tuple[HomedRayReturn, ...]:
    if not returns:
        return ()
    range_centers = np.asarray(range_centers_km, dtype=float)
    direct_range_km = min(float(solution.group_range_km) for solution in returns if math.isfinite(solution.group_range_km))
    range_step_km = float(np.median(np.diff(range_centers))) if range_centers.size > 1 else 25.0
    threshold_km = direct_range_km + max(75.0, 2.0 * range_step_km)
    filtered = tuple(
        solution for solution in returns
        if math.isfinite(solution.group_range_km) and float(solution.group_range_km) >= threshold_km
    )
    return filtered if filtered else tuple(returns)


def _ray_profile_along_track_km(case: TopsideCase, ray: object) -> tuple[np.ndarray, np.ndarray]:
    tx_ground = GeoPoint(case.tx_points[1].lat_deg, case.tx_points[1].lon_deg, 0.0)
    rx = case.rx_points[1]
    rx_ground = GeoPoint(rx.lat_deg, rx.lon_deg, 0.0)
    track_bearing_deg = initial_bearing_deg(tx_ground, rx_ground)
    lats = np.asarray(ray.path.get("lat", []), dtype=float)
    lons = np.asarray(ray.path.get("lon", []), dtype=float)
    heights = np.asarray(ray.path.get("height", []), dtype=float)
    valid = np.isfinite(lats) & np.isfinite(lons) & np.isfinite(heights) & (heights < 1e40)
    if np.count_nonzero(valid) < 2:
        return np.array([], dtype=float), np.array([], dtype=float)
    lats = lats[valid]
    lons = lons[valid]
    heights = heights[valid]
    ray_xyz = llh_to_ecef(lats, lons, heights * 1000.0)
    point_xyz = np.asarray(llh_to_ecef(rx.lat_deg, rx.lon_deg, rx.alt_km * 1000.0), dtype=float)
    if point_xyz.ndim > 1:
        point_xyz = point_xyz[0]
    start_index = 1
    if case.tx_points[1].alt_km > 400.0 and rx.alt_km > 400.0:
        threshold_alt_km = min(case.tx_points[1].alt_km, rx.alt_km) - 20.0
        below = np.flatnonzero(heights < threshold_alt_km)
        if below.size != 0:
            start_index = int(below[0])
    if start_index >= len(heights) - 1:
        return np.array([], dtype=float), np.array([], dtype=float)
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
        return np.array([], dtype=float), np.array([], dtype=float)
    best_index = int(np.argmin(distances))
    clip_index = start_index + best_index
    interp_lat = float(lats[clip_index] + t[best_index] * (lats[clip_index + 1] - lats[clip_index]))
    interp_lon = float(lons[clip_index] + t[best_index] * (lons[clip_index + 1] - lons[clip_index]))
    interp_height = float(heights[clip_index] + t[best_index] * (heights[clip_index + 1] - heights[clip_index]))
    clipped_lats = np.concatenate((lats[: clip_index + 1], np.array([interp_lat], dtype=float)))
    clipped_lons = np.concatenate((lons[: clip_index + 1], np.array([interp_lon], dtype=float)))
    clipped_heights = np.concatenate((heights[: clip_index + 1], np.array([interp_height], dtype=float)))
    along_km = np.asarray(_along_track_coordinate_km(tx_ground, track_bearing_deg, clipped_lats, clipped_lons), dtype=float)
    return np.maximum.accumulate(along_km), clipped_heights


def _profile_returns_in_window(
    returns: Sequence[HomedRayReturn],
    *,
    range_centers_km: np.ndarray,
    frequency_limits_mhz: tuple[float, float],
    max_miss_m: float,
) -> tuple[HomedRayReturn, ...]:
    filtered = _filter_ionospheric_returns(returns, range_centers_km)
    if not filtered:
        return ()
    fmin, fmax = sorted((float(frequency_limits_mhz[0]), float(frequency_limits_mhz[1])))
    return tuple(
        solution
        for solution in filtered
        if fmin - 1e-9 <= float(solution.frequency_mhz) <= fmax + 1e-9
        and math.isfinite(solution.miss_m)
        and float(solution.miss_m) <= float(max_miss_m)
    )


def _profile_curves_in_window(
    curves: Sequence[RayProfileCurve],
    *,
    range_centers_km: np.ndarray,
    frequency_limits_mhz: tuple[float, float],
    max_miss_m: float,
) -> tuple[RayProfileCurve, ...]:
    if not curves:
        return ()
    range_centers = np.asarray(range_centers_km, dtype=float)
    direct_range_km = min(float(curve.group_range_km) for curve in curves if math.isfinite(curve.group_range_km))
    range_step_km = float(np.median(np.diff(range_centers))) if range_centers.size > 1 else 25.0
    threshold_km = direct_range_km + max(75.0, 2.0 * range_step_km)
    fmin, fmax = sorted((float(frequency_limits_mhz[0]), float(frequency_limits_mhz[1])))
    filtered = tuple(
        curve
        for curve in curves
        if math.isfinite(curve.group_range_km)
        and float(curve.group_range_km) >= threshold_km
        and fmin - 1e-9 <= float(curve.frequency_mhz) <= fmax + 1e-9
        and math.isfinite(curve.miss_m)
        and float(curve.miss_m) <= float(max_miss_m)
    )
    return filtered if filtered else tuple(
        curve
        for curve in curves
        if fmin - 1e-9 <= float(curve.frequency_mhz) <= fmax + 1e-9
        and math.isfinite(curve.miss_m)
        and float(curve.miss_m) <= float(max_miss_m)
    )


def _render_frequency_plot_slice(
    problem: InverseProblem,
    case: TopsideCase,
    grid: IonosphereGrid,
    frequency_mhz: float,
) -> FrequencyPlotSlice:
    tracer = PointToPointRayTracer()
    center_tx = case.tx_points[1]
    center_rx = case.rx_points[1]
    nrange = problem.range_centers_km.size
    o_image = np.zeros((1, nrange), dtype=float)
    x_image = np.zeros((1, nrange), dtype=float)
    doppler_num = np.zeros((1, nrange), dtype=float)
    doppler_den = np.zeros((1, nrange), dtype=float)
    freq_axis = np.array([float(frequency_mhz)], dtype=float)
    o_curves: list[RayProfileCurve] = []
    x_curves: list[RayProfileCurve] = []
    for ox_mode in (1, -1):
        returns = _home_frequency_returns(
            tracer,
            tx=center_tx,
            rx=center_rx,
            grid=grid,
            fan_elevations_deg=case.fan_elevations_deg,
            fan_bearings_deg=case.fan_bearings_deg,
            frequency_mhz=float(frequency_mhz),
            ox_mode=ox_mode,
            config=problem.config,
        )
        target_image = o_image if ox_mode == 1 else x_image
        _accumulate_image(
            target_image,
            returns=returns,
            frequencies_mhz=freq_axis,
            range_edges_km=problem.range_edges_km,
        )
        returns_with_doppler = [
            replace(
                solution,
                doppler_hz=_estimate_return_doppler_hz(
                    tracer,
                    case=case,
                    center_return=solution,
                    grid=grid,
                    config=problem.config,
                ),
            )
            for solution in returns
        ]
        _accumulate_doppler(
            doppler_num,
            doppler_den,
            returns=returns_with_doppler,
            frequencies_mhz=freq_axis,
            range_edges_km=problem.range_edges_km,
        )
        target_curves = o_curves if ox_mode == 1 else x_curves
        for solution in returns:
            along_km, alt_km = _ray_profile_along_track_km(case, solution.ray)
            if along_km.size == 0 or alt_km.size == 0:
                continue
            target_curves.append(
                RayProfileCurve(
                    frequency_mhz=float(solution.frequency_mhz),
                    ox_mode=int(solution.ox_mode),
                    miss_m=float(solution.miss_m),
                    group_range_km=float(solution.group_range_km),
                    along_km=along_km,
                    alt_km=alt_km,
                )
            )
    del tracer
    gc.collect()
    return FrequencyPlotSlice(
        frequency_mhz=float(frequency_mhz),
        o_image=o_image[0],
        x_image=x_image[0],
        doppler_num=doppler_num[0],
        doppler_den=doppler_den[0],
        o_curves=tuple(o_curves),
        x_curves=tuple(x_curves),
    )


def _render_case_plot_data_chunked(
    problem: InverseProblem,
    case: TopsideCase,
    background_grid: IonosphereGrid,
    params: IonosphereFitParams,
) -> tuple[IonosphereGrid, CaseObservables, tuple[RayProfileCurve, ...], tuple[RayProfileCurve, ...]]:
    grid = _apply_fit_params_to_grid(problem, background_grid, params)
    total_shape = (problem.frequencies_mhz.size, problem.range_centers_km.size)
    o_image = np.zeros(total_shape, dtype=float)
    x_image = np.zeros(total_shape, dtype=float)
    doppler_num = np.zeros(total_shape, dtype=float)
    doppler_den = np.zeros(total_shape, dtype=float)
    o_curves: list[RayProfileCurve] = []
    x_curves: list[RayProfileCurve] = []
    for freq_index, frequency_mhz in enumerate(problem.frequencies_mhz):
        result = _render_frequency_plot_slice(problem, case, grid, float(frequency_mhz))
        o_image[freq_index] = result.o_image
        x_image[freq_index] = result.x_image
        doppler_num[freq_index] = result.doppler_num
        doppler_den[freq_index] = result.doppler_den
        o_curves.extend(result.o_curves)
        x_curves.extend(result.x_curves)
        del result
        gc.collect()

    sigma = (problem.config.blur_sigma_frequency_bins, problem.config.blur_sigma_range_bins)
    o_image = gaussian_filter(o_image, sigma=sigma, mode="nearest")
    x_image = gaussian_filter(x_image, sigma=sigma, mode="nearest")
    doppler_num = gaussian_filter(doppler_num, sigma=sigma, mode="nearest")
    doppler_den = gaussian_filter(doppler_den, sigma=sigma, mode="nearest")
    total_image = o_image + x_image
    total_scale = max(float(np.max(total_image)), 1e-9)
    o_norm = o_image / total_scale
    x_norm = x_image / total_scale
    total_norm = total_image / total_scale
    ox_split = (o_image - x_image) / np.maximum(total_image, 1e-6)
    doppler_image_hz = np.divide(doppler_num, doppler_den, out=np.zeros_like(doppler_num), where=doppler_den > 1e-9)
    observables = CaseObservables(
        name=case.name,
        kind=case.kind,
        total_image=total_norm,
        o_image=o_norm,
        x_image=x_norm,
        ox_split_image=ox_split,
        doppler_image_hz=doppler_image_hz,
        doppler_weight=np.clip(total_norm, 0.0, 1.0),
    )
    return grid, observables, tuple(o_curves), tuple(x_curves)


def plot_fit_summary(
    path: Path,
    demo: InverseDemoResult,
    *,
    case_name: str | None = None,
    plot_frequencies_mhz: np.ndarray | None = None,
) -> Path:
    import matplotlib.pyplot as plt  # type: ignore
    from matplotlib.cm import ScalarMappable  # type: ignore
    from matplotlib.colors import TwoSlopeNorm  # type: ignore

    case_index = _selected_case_index(demo, case_name=case_name)
    case = demo.problem.cases[case_index]
    background_grid = demo.problem.background_grids[case_index]
    if plot_frequencies_mhz is None:
        observed = demo.observed.cases[case_index]
        fitted = demo.fitted.cases[case_index]
        frequencies = np.asarray(demo.observed.frequencies_mhz, dtype=float)
        current_problem = demo.problem
        current_range_centers_km = np.asarray(demo.observed.range_centers_km, dtype=float)
        current_range_edges_km = np.asarray(demo.observed.range_edges_km, dtype=float)
        truth_grid, observed_o_curves, observed_x_curves = _trace_case_profile_curves(
            demo.problem,
            case,
            background_grid,
            demo.fit.truth_params,
        )
        fitted_grid, fitted_o_curves, fitted_x_curves = _trace_case_profile_curves(
            demo.problem,
            case,
            background_grid,
            demo.fit.fitted_params,
        )
        xmin, xmax = _dense_plot_xlim_mhz(
            frequencies,
            demo.observed.range_centers_km,
            observed.total_image,
            fitted.total_image,
        )
    else:
        frequencies = np.asarray(plot_frequencies_mhz, dtype=float)
        current_problem = _case_problem_with_frequencies(demo.problem, frequencies, range_bin_km=3.0)
        current_range_centers_km = np.asarray(current_problem.range_centers_km, dtype=float)
        current_range_edges_km = np.asarray(current_problem.range_edges_km, dtype=float)
        truth_grid, observed, observed_o_curves, observed_x_curves = _render_case_plot_data_chunked(
            current_problem,
            case,
            background_grid,
            demo.fit.truth_params,
        )
        fitted_grid, fitted, fitted_o_curves, fitted_x_curves = _render_case_plot_data_chunked(
            current_problem,
            case,
            background_grid,
            demo.fit.fitted_params,
        )
        xstep = float(np.median(np.diff(frequencies))) if frequencies.size > 1 else 0.5
        xmin = float(frequencies[0] - 0.5 * xstep)
        xmax = float(frequencies[-1] + 0.5 * xstep)
    observed_fmin, observed_fmax, observed_rmin, observed_rmax = _observed_support_extents(
        observed.total_image,
        frequencies,
        current_range_edges_km,
    )
    xmin = max(xmin, observed_fmin - 1.0)
    xmax = min(xmax, observed_fmax + 1.0)
    upper_range_min_km = max(float(current_range_edges_km[0]), observed_rmin - 100.0)
    upper_range_max_km = min(float(current_range_edges_km[-1]), observed_rmax + 100.0)
    extent = (
        float(frequencies[0] - 0.5 * (np.median(np.diff(frequencies)) if frequencies.size > 1 else 0.5)),
        float(frequencies[-1] + 0.5 * (np.median(np.diff(frequencies)) if frequencies.size > 1 else 0.5)),
        float(current_range_edges_km[-1]),
        float(current_range_edges_km[0]),
    )
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 8.9), constrained_layout=True)

    doppler_support = np.maximum(np.asarray(observed.doppler_weight, dtype=float), np.asarray(fitted.doppler_weight, dtype=float))
    _, observed_visible = _masked_visible_field(observed.doppler_image_hz, doppler_support)
    _, fitted_visible = _masked_visible_field(fitted.doppler_image_hz, doppler_support)
    visible_doppler = (
        np.concatenate([observed_visible, fitted_visible])
        if (observed_visible.size or fitted_visible.size)
        else np.array([0.0], dtype=float)
    )
    doppler_limit = max(float(np.nanmax(np.abs(visible_doppler))), 1.0)
    doppler_cmap = plt.get_cmap("RdBu_r")
    doppler_norm = TwoSlopeNorm(vcenter=0.0, vmin=-doppler_limit, vmax=doppler_limit)

    def _doppler_rgba(field: np.ndarray, support: np.ndarray) -> np.ndarray:
        rgba = doppler_cmap(doppler_norm(np.asarray(field, dtype=float).T))
        alpha = np.clip(np.asarray(support, dtype=float).T, 0.0, 1.0)
        alpha = np.where(np.isfinite(np.asarray(field, dtype=float).T), alpha, 0.0)
        rgba[..., 3] = alpha
        return rgba

    ionogram_fields = (
        (axes[0, 0], observed.doppler_image_hz, observed.total_image, "Observed"),
        (axes[0, 1], fitted.doppler_image_hz, fitted.total_image, "Fitted"),
    )
    for ax, field, support, title in ionogram_fields:
        ax.imshow(
            _doppler_rgba(field, support),
            origin="upper",
            aspect="auto",
            extent=extent,
            interpolation="nearest",
        )
        ax.set_title(title)

    tx_ground, rx_ground, track_lats_deg, track_lons_deg, along_track_km = _case_track(case)
    total_ground_track_km = float(along_track_km[-1])
    swath_truth_pf_mhz = _profile_swath_plasma_frequency_mhz(truth_grid, track_lats_deg, track_lons_deg)
    swath_fitted_pf_mhz = _profile_swath_plasma_frequency_mhz(fitted_grid, track_lats_deg, track_lons_deg)
    altitudes_km = np.asarray(truth_grid.altitudes_km, dtype=float)
    plotted_altitudes_km, plotted_truth_swath = _resample_vertical_display(
        swath_truth_pf_mhz,
        altitudes_km,
        step_km=3.0,
    )
    _, plotted_fitted_swath = _resample_vertical_display(
        swath_fitted_pf_mhz,
        altitudes_km,
        step_km=3.0,
    )
    along_edges = _centers_to_edges(along_track_km)
    along_edges[0] = -50.0
    along_edges[-1] = total_ground_track_km + 50.0
    alt_edges = _centers_to_edges(plotted_altitudes_km)
    plotted_alt_edges = alt_edges
    if alt_edges.size >= 2 and alt_edges[0] > 0.0:
        zero_row = np.zeros((1, plotted_truth_swath.shape[1]), dtype=float)
        plotted_truth_swath = np.vstack((zero_row, plotted_truth_swath))
        plotted_fitted_swath = np.vstack((zero_row, plotted_fitted_swath))
        plotted_alt_edges = np.concatenate(([0.0], alt_edges))
    along_mesh, alt_mesh = np.meshgrid(along_edges, plotted_alt_edges, indexing="xy")
    visible_profile = np.concatenate((plotted_truth_swath.ravel(), plotted_fitted_swath.ravel()))
    visible_profile = visible_profile[np.isfinite(visible_profile)]
    pf_vmax = float(np.max(visible_profile)) if visible_profile.size else 1.0
    if pf_vmax <= 0.0:
        pf_vmax = 1.0

    ordinary_color = "white"
    extraordinary_color = "white"
    profile_fields = (
        (
            axes[1, 0],
            plotted_truth_swath,
            _profile_curves_in_window(
                observed_o_curves,
                range_centers_km=current_range_centers_km,
                frequency_limits_mhz=(xmin, xmax),
                max_miss_m=demo.problem.config.homing_tolerance_m,
            ),
            _profile_curves_in_window(
                observed_x_curves,
                range_centers_km=current_range_centers_km,
                frequency_limits_mhz=(xmin, xmax),
                max_miss_m=demo.problem.config.homing_tolerance_m,
            ),
            "",
        ),
        (
            axes[1, 1],
            plotted_fitted_swath,
            _profile_curves_in_window(
                fitted_o_curves,
                range_centers_km=current_range_centers_km,
                frequency_limits_mhz=(xmin, xmax),
                max_miss_m=demo.problem.config.homing_tolerance_m,
            ),
            _profile_curves_in_window(
                fitted_x_curves,
                range_centers_km=current_range_centers_km,
                frequency_limits_mhz=(xmin, xmax),
                max_miss_m=demo.problem.config.homing_tolerance_m,
            ),
            "",
        ),
    )
    profile_images = []
    for ax, swath_pf_mhz, ordinary_returns, extraordinary_returns, title in profile_fields:
        image = ax.pcolormesh(
            along_mesh,
            alt_mesh,
            swath_pf_mhz,
            cmap="viridis",
            vmin=0.0,
            vmax=pf_vmax,
            shading="flat",
            antialiased=False,
            rasterized=True,
        )
        profile_images.append(image)
        for curve in ordinary_returns:
            if curve.along_km.size:
                ax.plot(curve.along_km, curve.alt_km, color=ordinary_color, linewidth=1.25, alpha=0.9)
        for curve in extraordinary_returns:
            if curve.along_km.size:
                ax.plot(curve.along_km, curve.alt_km, color=extraordinary_color, linewidth=1.1, alpha=0.9, linestyle="--")
        ax.scatter(
            [0.0, total_ground_track_km],
            [case.tx_points[1].alt_km, case.rx_points[1].alt_km],
            c=["gold", "tab:red"],
            edgecolors="black",
            linewidths=0.7,
            s=42,
            zorder=5,
        )
        ax.set_title(title)

    for row_index, row in enumerate(axes):
        for col_index, ax in enumerate(row):
            ax.set_facecolor("black" if row_index == 0 else "white")
            ax.set_axisbelow(True)
            ax.minorticks_on()
            if row_index == 0:
                ax.grid(True, which="major", color="white", alpha=0.12, linewidth=0.7)
                ax.grid(True, which="minor", color="white", alpha=0.06, linewidth=0.5, linestyle=":")
            else:
                ax.grid(True, which="major", alpha=0.22, linewidth=0.7)
                ax.grid(True, which="minor", alpha=0.10, linewidth=0.5, linestyle=":")
            if row_index == 0:
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(upper_range_max_km, upper_range_min_km)
                ax.set_xlabel("Frequency (MHz)")
                if col_index == 0:
                    ax.set_ylabel("Virtual range (km)")
            else:
                ax.set_xlim(-50.0, total_ground_track_km + 50.0)
                ax.set_ylim(0.0, max(float(plotted_altitudes_km[-1]), case.tx_points[1].alt_km, case.rx_points[1].alt_km) + 10.0)
                ax.set_xlabel("Ground-track distance (km)")
                if col_index == 0:
                    ax.set_ylabel("Altitude (km)")
    axes[0, 0].tick_params(labelbottom=True, color="white", labelcolor="black")
    axes[0, 1].tick_params(labelbottom=True, color="white", labelcolor="black")
    axes[0, 0].tick_params(labelleft=True, color="white", labelcolor="black")
    axes[0, 1].tick_params(labelleft=False, color="white", labelcolor="black")
    axes[0, 0].xaxis.label.set_color("black")
    axes[0, 1].xaxis.label.set_color("black")
    axes[0, 0].yaxis.label.set_color("black")
    axes[0, 0].title.set_color("black")
    axes[0, 1].title.set_color("black")
    for ax in axes[0, :]:
        for spine in ax.spines.values():
            spine.set_color("white")

    doppler_cbar = fig.colorbar(
        ScalarMappable(norm=doppler_norm, cmap=doppler_cmap),
        ax=axes[0, :],
        pad=0.02,
        shrink=0.92,
    )
    doppler_cbar.set_label("Doppler (Hz); opacity set by normalized power")
    profile_cbar = fig.colorbar(profile_images[0], ax=axes[1, :], pad=0.02, shrink=0.92)
    profile_cbar.set_label("Plasma frequency (MHz)")
    axes[1, 1].legend(
        handles=[
            plt.Line2D([0.0], [0.0], color=ordinary_color, linewidth=1.4, label="O-mode"),
            plt.Line2D([0.0], [0.0], color=extraordinary_color, linewidth=1.2, linestyle="--", label="X-mode"),
        ],
        loc="upper right",
        framealpha=0.9,
    )

    fig.suptitle(
        f"{demo.problem.cases[case_index].name}: truth={demo.fit.truth_cost:.4f}, fitted={demo.fit.fitted_cost:.4f}",
        fontsize=14,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


def _case_track(case: TopsideCase) -> tuple[GeoPoint, GeoPoint, np.ndarray, np.ndarray, np.ndarray]:
    tx_ground = GeoPoint(case.tx_points[1].lat_deg, case.tx_points[1].lon_deg, 0.0)
    rx_ground = GeoPoint(case.rx_points[1].lat_deg, case.rx_points[1].lon_deg, 0.0)
    if case.kind == "vertical":
        track_lats = np.array([tx_ground.lat_deg, tx_ground.lat_deg + 0.25], dtype=float)
        track_lons = np.array([tx_ground.lon_deg, tx_ground.lon_deg], dtype=float)
        along_track_km = np.array([0.0, 1.0], dtype=float)
    else:
        track_lats, track_lons = great_circle_waypoints(tx_ground, rx_ground, count=256)
        along_track_km = np.linspace(0.0, _surface_distance_km(tx_ground, rx_ground), len(track_lats), dtype=float)
    return tx_ground, rx_ground, np.asarray(track_lats, dtype=float), np.asarray(track_lons, dtype=float), along_track_km


def _profile_swath_plasma_frequency_mhz(
    grid: IonosphereGrid,
    track_lats_deg: np.ndarray,
    track_lons_deg: np.ndarray,
) -> np.ndarray:
    density_m3 = np.asarray(grid.iono_en_grid, dtype=float) * 1e6
    sample_lons_deg = np.asarray(
        [coerce_longitude_for_grid(float(lon), grid.longitudes_deg) for lon in np.asarray(track_lons_deg, dtype=float)],
        dtype=float,
    )
    lat_mesh, alt_mesh = np.meshgrid(np.asarray(track_lats_deg, dtype=float), np.asarray(grid.altitudes_km, dtype=float), indexing="xy")
    lon_mesh, _ = np.meshgrid(sample_lons_deg, np.asarray(grid.altitudes_km, dtype=float), indexing="xy")
    points = np.column_stack((lat_mesh.ravel(), lon_mesh.ravel(), alt_mesh.ravel()))
    interpolator = RegularGridInterpolator(
        (grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km),
        density_m3,
        bounds_error=False,
        fill_value=np.nan,
    )
    density_swath_m3 = interpolator(points).reshape(len(grid.altitudes_km), len(track_lats_deg))
    density_swath_m3 = np.where(np.isfinite(density_swath_m3) & (density_swath_m3 > 0.0), density_swath_m3, 0.0)
    return 8.98e-6 * np.sqrt(np.maximum(density_swath_m3, 0.0))


def plot_ionosphere_comparison(
    path: Path,
    demo: InverseDemoResult,
    *,
    case_name: str | None = None,
) -> Path:
    import matplotlib.pyplot as plt  # type: ignore

    if demo.problem.global_background_grid is None:
        raise ValueError("global background grid was discarded; rerun with keep_global_background=True")
    case_index = _selected_case_index(demo, case_name=case_name)
    case = demo.problem.cases[case_index]
    truth_grid = _apply_fit_params_to_grid(demo.problem, demo.problem.global_background_grid, demo.fit.truth_params)
    fitted_grid = _apply_fit_params_to_grid(demo.problem, demo.problem.global_background_grid, demo.fit.fitted_params)

    tx_ground, rx_ground, track_lats_deg, track_lons_deg, along_track_km = _case_track(case)
    plot_track_lons_deg = np.asarray(
        [coerce_longitude_for_grid(float(lon), truth_grid.longitudes_deg) for lon in track_lons_deg],
        dtype=float,
    )
    tx_lon_plot = coerce_longitude_for_grid(tx_ground.lon_deg, truth_grid.longitudes_deg)
    rx_lon_plot = coerce_longitude_for_grid(rx_ground.lon_deg, truth_grid.longitudes_deg)

    truth_peak_density_m3 = np.max(np.asarray(truth_grid.iono_en_grid, dtype=float) * 1e6, axis=2)
    fitted_peak_density_m3 = np.max(np.asarray(fitted_grid.iono_en_grid, dtype=float) * 1e6, axis=2)
    truth_fof2_mhz = 8.98e-6 * np.sqrt(np.maximum(truth_peak_density_m3, 0.0))
    fitted_fof2_mhz = 8.98e-6 * np.sqrt(np.maximum(fitted_peak_density_m3, 0.0))

    lon_edges = _centers_to_edges(np.asarray(truth_grid.longitudes_deg, dtype=float))
    lat_edges = _centers_to_edges(np.asarray(truth_grid.latitudes_deg, dtype=float))
    lon_mesh, lat_mesh = np.meshgrid(lon_edges, lat_edges, indexing="xy")

    lat_values = [tx_ground.lat_deg, rx_ground.lat_deg, *track_lats_deg.tolist()]
    lon_values = [tx_lon_plot, rx_lon_plot, *plot_track_lons_deg.tolist()]
    lat_min = min(lat_values)
    lat_max = max(lat_values)
    lon_min = min(lon_values)
    lon_max = max(lon_values)
    lat_margin = max(1.0, 0.08 * max(lat_max - lat_min, 0.5))
    lon_margin = max(1.0, 0.08 * max(lon_max - lon_min, 0.5))
    xlim_plan = (lon_min - lon_margin, lon_max + lon_margin)
    ylim_plan = (lat_min - lat_margin, lat_max + lat_margin)

    combined_plan = np.maximum(truth_fof2_mhz, fitted_fof2_mhz)
    plan_vmin, plan_vmax = _visible_color_limits(
        combined_plan,
        x_centers=np.asarray(truth_grid.longitudes_deg, dtype=float),
        y_centers=np.asarray(truth_grid.latitudes_deg, dtype=float),
        x_limits=xlim_plan,
        y_limits=ylim_plan,
    )

    truth_profile_pf_mhz = _profile_swath_plasma_frequency_mhz(truth_grid, track_lats_deg, track_lons_deg)
    fitted_profile_pf_mhz = _profile_swath_plasma_frequency_mhz(fitted_grid, track_lats_deg, track_lons_deg)
    altitudes_km = np.asarray(truth_grid.altitudes_km, dtype=float)
    along_edges = _centers_to_edges(along_track_km)
    alt_edges = _centers_to_edges(altitudes_km)
    profile_ymin = 0.0
    profile_ymax = max(float(altitudes_km[-1]), case.tx_points[1].alt_km, case.rx_points[1].alt_km)
    profile_ymax += 10.0
    plotted_truth_profile_pf_mhz = truth_profile_pf_mhz
    plotted_fitted_profile_pf_mhz = fitted_profile_pf_mhz
    plotted_alt_edges = alt_edges
    if alt_edges.size >= 2 and alt_edges[0] > 0.0:
        zero_row = np.zeros((1, truth_profile_pf_mhz.shape[1]), dtype=float)
        plotted_truth_profile_pf_mhz = np.vstack((zero_row, truth_profile_pf_mhz))
        plotted_fitted_profile_pf_mhz = np.vstack((zero_row, fitted_profile_pf_mhz))
        plotted_alt_edges = np.concatenate(([0.0], alt_edges))
    along_mesh, alt_mesh = np.meshgrid(along_edges, plotted_alt_edges, indexing="xy")

    row_visible = (plotted_alt_edges[:-1] < profile_ymax) & (plotted_alt_edges[1:] > profile_ymin)
    col_visible = (along_edges[:-1] < along_track_km[-1]) & (along_edges[1:] > along_track_km[0])
    visible_truth = plotted_truth_profile_pf_mhz[np.ix_(row_visible, col_visible)]
    visible_fitted = plotted_fitted_profile_pf_mhz[np.ix_(row_visible, col_visible)]
    visible_profile = np.concatenate((visible_truth.ravel(), visible_fitted.ravel()))
    visible_profile = visible_profile[np.isfinite(visible_profile)]
    if visible_profile.size:
        profile_vmin = 0.0
        profile_vmax = float(np.max(visible_profile))
        if profile_vmax <= profile_vmin:
            profile_vmax = profile_vmin + 1e-6
    else:
        profile_vmin = 0.0
        profile_vmax = 1.0

    fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.0))
    plan_fields = (
        (axes[0, 0], truth_fof2_mhz, "Truth foF2 / peak plasma frequency"),
        (axes[0, 1], fitted_fof2_mhz, "Retrieved foF2 / peak plasma frequency"),
    )
    for ax, field, title in plan_fields:
        surface = ax.pcolormesh(
            lon_mesh,
            lat_mesh,
            field,
            cmap="viridis",
            vmin=plan_vmin,
            vmax=plan_vmax,
            shading="flat",
            antialiased=False,
            rasterized=True,
        )
        ax.plot(plot_track_lons_deg, track_lats_deg, color="white", linewidth=2.2, alpha=0.95)
        ax.scatter([tx_lon_plot], [tx_ground.lat_deg], color="gold", edgecolors="black", linewidths=0.8, s=120, marker="*", zorder=5)
        ax.scatter([rx_lon_plot], [rx_ground.lat_deg], color="tab:red", edgecolors="white", linewidths=0.8, s=52, zorder=5)
        ax.set_xlim(*xlim_plan)
        ax.set_ylim(*ylim_plan)
        ax.set_xlabel("Longitude (deg)")
        ax.set_ylabel("Latitude (deg)")
        ax.set_title(title)
        fig.colorbar(surface, ax=ax, pad=0.02, label="foF2 / peak plasma frequency (MHz)")

    profile_fields = (
        (axes[1, 0], plotted_truth_profile_pf_mhz, "Truth altitude profile"),
        (axes[1, 1], plotted_fitted_profile_pf_mhz, "Retrieved altitude profile"),
    )
    total_ground_track_km = float(along_track_km[-1])
    for ax, field, title in profile_fields:
        surface = ax.pcolormesh(
            along_mesh,
            alt_mesh,
            field,
            cmap="viridis",
            vmin=profile_vmin,
            vmax=profile_vmax,
            shading="flat",
            antialiased=False,
            rasterized=True,
        )
        ax.scatter(
            [0.0, total_ground_track_km],
            [case.tx_points[1].alt_km, case.rx_points[1].alt_km],
            c=["gold", "tab:red"],
            edgecolors="black",
            linewidths=0.8,
            s=54,
            zorder=5,
        )
        ax.set_xlim(0.0, max(total_ground_track_km, 1.0))
        ax.set_ylim(profile_ymin, profile_ymax)
        ax.set_xlabel("Ground-track distance (km)")
        ax.set_ylabel("Altitude (km)")
        ax.set_title(title)
        fig.colorbar(surface, ax=ax, pad=0.02, label="Plasma frequency (MHz)")

    fig.suptitle(
        f"{case.name}: truth vs retrieved ionosphere",
        fontsize=14,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    return path


def build_parser() -> argparse.ArgumentParser:
    config = TopsideInverseConfig()
    parser = argparse.ArgumentParser(
        description="Synthetic multisatellite topside inverse test case using a fixed PHaRLAP launch fan and AMPERE-derived spoof orbits.",
    )
    parser.add_argument("--ampere-file", type=Path, default=config.ampere_file)
    parser.add_argument("--center-time", default=config.center_utc.isoformat())
    parser.add_argument("--planes", nargs="+", type=int, default=list(config.planes))
    parser.add_argument("--spoof-alt-km", type=float, default=config.spoof_altitude_km)
    parser.add_argument("--oblique-separation-km", type=float, default=config.oblique_separation_km)
    parser.add_argument("--frequencies", nargs="+", type=float, default=list(config.frequencies_mhz))
    parser.add_argument("--solver-maxiter", type=int, default=config.solver_maxiter)
    parser.add_argument("--solver-popsize", type=int, default=config.solver_popsize)
    parser.add_argument("--solver-seed", type=int, default=config.solver_seed)
    parser.add_argument(
        "--grid-cache",
        type=Path,
        default=None,
        help="Optional netCDF file or directory for the cached global background grid.",
    )
    parser.add_argument("--rebuild-grid", action="store_true")
    parser.add_argument("--plot-out", type=Path, default=None)
    parser.add_argument("--ionosphere-plot-out", type=Path, default=None)
    parser.add_argument("--fit-out", type=Path, default=None)
    parser.add_argument("--case-name", default=None)
    parser.add_argument(
        "--fit-in",
        type=Path,
        default=None,
        help="Load previously saved fit/demo JSON and skip the inverse solve.",
    )
    parser.add_argument(
        "--fit-case-name",
        action="append",
        default=None,
        help="Restrict the inverse solve to one or more named cases before plotting.",
    )
    parser.add_argument("--plot-frequency-start-mhz", type=float, default=None)
    parser.add_argument("--plot-frequency-stop-mhz", type=float, default=None)
    parser.add_argument("--plot-frequency-step-khz", type=float, default=None)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    config = replace(
        TopsideInverseConfig(),
        ampere_file=args.ampere_file.expanduser(),
        center_utc=_parse_time(args.center_time),
        planes=tuple(int(value) for value in args.planes),
        spoof_altitude_km=float(args.spoof_alt_km),
        oblique_separation_km=float(args.oblique_separation_km),
        frequencies_mhz=tuple(float(value) for value in args.frequencies),
        solver_maxiter=int(args.solver_maxiter),
        solver_popsize=int(args.solver_popsize),
        solver_seed=int(args.solver_seed),
        grid_cache_path=None if args.grid_cache is None else args.grid_cache.expanduser(),
        rebuild_grid=bool(args.rebuild_grid),
    )
    try:
        if args.fit_in is None and args.fit_out is not None and args.plot_out is None and args.ionosphere_plot_out is None:
            fit = run_inverse_fit_only(
                config,
                fit_case_names=args.fit_case_name,
                keep_global_background=False,
            )
            payload = {"fit": fit.to_dict()}
            args.fit_out.expanduser().write_text(json.dumps(payload, indent=2))
            print(json.dumps(payload, indent=2))
            return
        if args.fit_in is None:
            demo = run_inverse_demo(
                config,
                fit_case_names=args.fit_case_name,
                keep_global_background=bool(args.ionosphere_plot_out),
            )
        else:
            fit = _load_fit_result_json(args.fit_in.expanduser())
            problem = _subset_problem_cases(build_inverse_problem(config), args.fit_case_name)
            if not args.ionosphere_plot_out:
                problem = replace(problem, global_background_grid=None)
                gc.collect()
            empty_cases = tuple(_empty_case_observables(problem, case) for case in problem.cases)
            empty_dataset = SyntheticDataset(
                frequencies_mhz=problem.frequencies_mhz,
                range_edges_km=problem.range_edges_km,
                range_centers_km=problem.range_centers_km,
                cases=empty_cases,
            )
            demo = InverseDemoResult(
                problem=problem,
                observed=empty_dataset,
                fitted=empty_dataset,
                fit=fit,
            )
    except PyLapImportError as exc:
        raise SystemExit(str(exc)) from exc
    plot_frequencies_mhz = None
    if args.plot_frequency_step_khz is not None:
        start_mhz = float(args.plot_frequency_start_mhz) if args.plot_frequency_start_mhz is not None else float(min(config.frequencies_mhz))
        stop_mhz = float(args.plot_frequency_stop_mhz) if args.plot_frequency_stop_mhz is not None else float(max(config.frequencies_mhz))
        plot_frequencies_mhz = make_frequency_axis_mhz(start_mhz, stop_mhz, float(args.plot_frequency_step_khz))
    if args.plot_out is not None:
        plot_fit_summary(
            args.plot_out.expanduser(),
            demo,
            case_name=args.case_name,
            plot_frequencies_mhz=plot_frequencies_mhz,
        )
    if args.ionosphere_plot_out is not None:
        plot_ionosphere_comparison(
            args.ionosphere_plot_out.expanduser(),
            demo,
            case_name=args.case_name,
        )
    print(json.dumps(demo.to_dict(), indent=2))


if __name__ == "__main__":
    main()
