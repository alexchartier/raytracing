from __future__ import annotations

import argparse
import datetime as dt
import json
import math
from dataclasses import dataclass

import numpy as np

if __package__:
    from .geometry import (
        GeoPoint,
        WGS84_A_M,
        WGS84_E2,
        initial_bearing_deg,
        llh_to_ecef,
        make_regional_lat_lon_grids,
        relaz_to_ecef_unit,
        wrap_longitude,
    )
    from .grid import IonosphereGrid, build_pyiri_grid
    from .tracer import PointToPointRayTracer, PyLapImportError, RayTrace
else:
    from geometry import (
        GeoPoint,
        WGS84_A_M,
        WGS84_E2,
        initial_bearing_deg,
        llh_to_ecef,
        make_regional_lat_lon_grids,
        relaz_to_ecef_unit,
        wrap_longitude,
    )
    from grid import IonosphereGrid, build_pyiri_grid
    from tracer import PointToPointRayTracer, PyLapImportError, RayTrace


@dataclass(frozen=True)
class GroundToSpaceScenario:
    when: dt.datetime
    tx: GeoPoint
    rx: GeoPoint
    frequency_mhz: float = 12.0
    f107: float = 120.0
    ap_daily: float = 8.0
    nhops: int = 1
    homing_tolerance_m: float = 500.0


@dataclass(frozen=True)
class GroundToSpaceDemoResult:
    mode: str
    scenario: GroundToSpaceScenario
    ray: RayTrace
    slant_range_km: float
    line_of_sight_elevation_deg: float
    line_of_sight_bearing_deg: float
    max_path_alt_km: float

    def to_dict(self) -> dict:
        return {
            "mode": self.mode,
            "when": self.scenario.when.isoformat(),
            "tx": {
                "lat_deg": self.scenario.tx.lat_deg,
                "lon_deg": self.scenario.tx.lon_deg,
                "alt_km": self.scenario.tx.alt_km,
            },
            "rx": {
                "lat_deg": self.scenario.rx.lat_deg,
                "lon_deg": self.scenario.rx.lon_deg,
                "alt_km": self.scenario.rx.alt_km,
            },
            "line_of_sight": {
                "slant_range_km": self.slant_range_km,
                "elevation_deg": self.line_of_sight_elevation_deg,
                "bearing_deg": self.line_of_sight_bearing_deg,
            },
            "solution": {
                **self.ray.to_dict(),
                "max_path_alt_km": self.max_path_alt_km,
            },
        }


def default_scenario() -> GroundToSpaceScenario:
    return GroundToSpaceScenario(
        when=dt.datetime(2020, 1, 15, 12, 0, 0),
        tx=GeoPoint(34.873, -106.614, 1.6),
        rx=GeoPoint(32.540, -102.450, 550.0),
    )


def _parse_time(text: str) -> dt.datetime:
    when = dt.datetime.fromisoformat(text)
    if when.tzinfo is not None:
        when = when.astimezone(dt.timezone.utc).replace(tzinfo=None)
    return when


def _grid_parms(latitudes_deg: np.ndarray, longitudes_deg: np.ndarray, altitudes_km: np.ndarray) -> list[float]:
    lat_step = float(latitudes_deg[1] - latitudes_deg[0]) if len(latitudes_deg) > 1 else 0.0
    lon_step = float(longitudes_deg[1] - longitudes_deg[0]) if len(longitudes_deg) > 1 else 0.0
    alt_step = float(altitudes_km[1] - altitudes_km[0]) if len(altitudes_km) > 1 else 0.0
    return [
        float(latitudes_deg[0]), lat_step, float(len(latitudes_deg)),
        float(longitudes_deg[0]), lon_step, float(len(longitudes_deg)),
        float(altitudes_km[0]), alt_step, float(len(altitudes_km)),
    ]


def _ecef_to_llh(points_xyz_m: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    points = np.asarray(points_xyz_m, dtype=float)
    x = points[..., 0]
    y = points[..., 1]
    z = points[..., 2]
    lon = np.arctan2(y, x)
    p = np.hypot(x, y)

    lat = np.arctan2(z, p * (1.0 - WGS84_E2))
    alt = np.zeros_like(lat)
    for _ in range(6):
        sin_lat = np.sin(lat)
        cos_lat = np.cos(lat)
        prime_vertical = WGS84_A_M / np.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
        safe_cos = np.where(np.abs(cos_lat) < 1e-12, np.sign(cos_lat) * 1e-12 + (cos_lat == 0.0) * 1e-12, cos_lat)
        alt = p / safe_cos - prime_vertical
        lat = np.arctan2(z, p * (1.0 - WGS84_E2 * prime_vertical / np.maximum(prime_vertical + alt, 1.0)))

    return np.degrees(lat), np.degrees(lon), alt


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


def build_synthetic_demo_grid(tx: GeoPoint, rx: GeoPoint) -> IonosphereGrid:
    latitudes_deg, longitudes_deg = make_regional_lat_lon_grids(
        tx,
        rx,
        lat_step_deg=2.0,
        lon_step_deg=2.0,
        lat_margin_deg=2.0,
        lon_margin_deg=2.0,
        waypoint_count=32,
    )
    if len(latitudes_deg) < 2:
        latitudes_deg = np.array([tx.lat_deg - 1.0, tx.lat_deg + 1.0], dtype=float)
    if len(longitudes_deg) < 2:
        longitudes_deg = np.array([tx.lon_deg - 1.0, tx.lon_deg + 1.0], dtype=float)

    altitudes_km = np.array([60.0, max(650.0, math.ceil(rx.alt_km / 50.0) * 50.0)], dtype=float)
    shape = (len(latitudes_deg), len(longitudes_deg), len(altitudes_km))
    zeros = np.zeros(shape, dtype=float)
    parms = _grid_parms(latitudes_deg, longitudes_deg, altitudes_km)
    return IonosphereGrid(
        latitudes_deg=latitudes_deg,
        longitudes_deg=longitudes_deg,
        altitudes_km=altitudes_km,
        iono_en_grid=zeros.copy(),
        iono_en_grid_5=zeros.copy(),
        collision_freq=zeros.copy(),
        iono_grid_parms=parms,
        Bx=zeros.copy(),
        By=zeros.copy(),
        Bz=zeros.copy(),
        geomag_grid_parms=parms.copy(),
        metadata={"mode": "synthetic_space_demo"},
    )


class StraightLineSpaceBackend:
    def __init__(self, target: GeoPoint, sample_count: int = 96, overshoot_factor: float = 1.25):
        self.target = target
        self.sample_count = sample_count
        self.overshoot_factor = overshoot_factor

    def trace(self, origin, elevations_deg, bearings_deg, freqs_mhz, ox_mode, nhops, tol, *, grid=None, state_vector=None):
        origin_xyz = np.asarray(llh_to_ecef(origin.lat_deg, origin.lon_deg, origin.alt_km * 1000.0), dtype=float)
        target_xyz = np.asarray(llh_to_ecef(self.target.lat_deg, self.target.lon_deg, self.target.alt_km * 1000.0), dtype=float)
        if origin_xyz.ndim > 1:
            origin_xyz = origin_xyz[0]
        if target_xyz.ndim > 1:
            target_xyz = target_xyz[0]

        target_slant_range_m = float(np.linalg.norm(target_xyz - origin_xyz))
        path_range_m = max(target_slant_range_m * self.overshoot_factor, target_slant_range_m + 150_000.0)
        samples = np.linspace(0.0, path_range_m, self.sample_count)

        rays: list[RayTrace] = []
        for elevation_deg, bearing_deg, frequency_mhz in zip(elevations_deg, bearings_deg, freqs_mhz):
            direction_xyz = np.asarray(
                relaz_to_ecef_unit(float(elevation_deg), float(bearing_deg), origin.lat_deg, origin.lon_deg),
                dtype=float,
            )
            if direction_xyz.ndim > 1:
                direction_xyz = direction_xyz[0]

            path_xyz = origin_xyz[None, :] + samples[:, None] * direction_xyz[None, :]
            lat_deg, lon_deg, alt_m = _ecef_to_llh(path_xyz)
            path = {
                "initial_elev": float(elevation_deg),
                "initial_bearing": float(bearing_deg),
                "frequency": float(frequency_mhz),
                "lat": np.asarray(lat_deg, dtype=float),
                "lon": np.asarray([wrap_longitude(value) for value in lon_deg], dtype=float),
                "height": np.asarray(alt_m / 1000.0, dtype=float),
                "group_range": np.asarray(samples / 1000.0, dtype=float),
                "geometric_distance": np.asarray(samples / 1000.0, dtype=float),
                "absorption": np.zeros(self.sample_count, dtype=float),
            }
            rays.append(
                RayTrace(
                    summary={"backend": "straight_line_space_demo", "nhops": int(nhops)},
                    path=path,
                    state={"target_slant_range_m": target_slant_range_m},
                )
            )
        return rays


def run_ground_to_space_demo(mode: str = "synthetic", scenario: GroundToSpaceScenario | None = None) -> GroundToSpaceDemoResult:
    scenario = scenario or default_scenario()
    los_bearing_deg, los_elevation_deg, slant_range_km = _line_of_sight_angles(scenario.tx, scenario.rx)

    if mode == "synthetic":
        grid = build_synthetic_demo_grid(scenario.tx, scenario.rx)
        tracer = PointToPointRayTracer(backend=StraightLineSpaceBackend(scenario.rx))
    elif mode == "real":
        grid = build_pyiri_grid(
            scenario.when,
            scenario.tx,
            scenario.rx,
            f107=scenario.f107,
            ap_daily=scenario.ap_daily,
            alt_min_km=60.0,
            alt_max_km=max(650.0, scenario.rx.alt_km + 100.0),
            alt_step_km=5.0,
            lat_step_deg=1.0,
            lon_step_deg=1.0,
            lat_margin_deg=5.0,
            lon_margin_deg=5.0,
            d_region_model="fpt2018",
        )
        tracer = PointToPointRayTracer()
    else:
        raise ValueError(f"unsupported mode: {mode}")

    ray = tracer.trace_link(
        tx=scenario.tx,
        rx=scenario.rx,
        frequency_mhz=scenario.frequency_mhz,
        grid=grid,
        nhops=scenario.nhops,
        homing_tolerance_m=scenario.homing_tolerance_m,
    )
    finite_heights = np.asarray(ray.path.get("height", []), dtype=float)
    finite_heights = finite_heights[np.isfinite(finite_heights) & (finite_heights < 1e40)]
    max_path_alt_km = float(np.max(finite_heights)) if finite_heights.size else math.nan
    return GroundToSpaceDemoResult(
        mode=mode,
        scenario=scenario,
        ray=ray,
        slant_range_km=slant_range_km,
        line_of_sight_elevation_deg=los_elevation_deg,
        line_of_sight_bearing_deg=los_bearing_deg,
        max_path_alt_km=max_path_alt_km,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Ground-to-space 3D homing example built on the existing python_raytrace solver.",
    )
    parser.add_argument("--mode", choices=("synthetic", "real"), default="synthetic")
    parser.add_argument("--time", default=default_scenario().when.isoformat(), help="UTC time in ISO-8601 format.")
    parser.add_argument("--tx", nargs=3, type=float, metavar=("LAT", "LON", "ALT_KM"), default=None)
    parser.add_argument("--rx", nargs=3, type=float, metavar=("LAT", "LON", "ALT_KM"), default=None)
    parser.add_argument("--freq", type=float, default=default_scenario().frequency_mhz, help="Frequency in MHz.")
    parser.add_argument("--f107", type=float, default=default_scenario().f107)
    parser.add_argument("--ap-daily", type=float, default=default_scenario().ap_daily)
    parser.add_argument("--nhops", type=int, default=default_scenario().nhops)
    parser.add_argument("--homing-tolerance-m", type=float, default=default_scenario().homing_tolerance_m)
    return parser


def _scenario_from_args(args: argparse.Namespace) -> GroundToSpaceScenario:
    defaults = default_scenario()
    tx = GeoPoint(*args.tx) if args.tx is not None else defaults.tx
    rx = GeoPoint(*args.rx) if args.rx is not None else defaults.rx
    return GroundToSpaceScenario(
        when=_parse_time(args.time),
        tx=tx,
        rx=rx,
        frequency_mhz=float(args.freq),
        f107=float(args.f107),
        ap_daily=float(args.ap_daily),
        nhops=int(args.nhops),
        homing_tolerance_m=float(args.homing_tolerance_m),
    )


def main() -> None:
    args = build_parser().parse_args()
    scenario = _scenario_from_args(args)
    try:
        result = run_ground_to_space_demo(args.mode, scenario=scenario)
    except PyLapImportError as exc:
        raise SystemExit(f"real mode requires a working pylap build: {exc}") from exc
    print(json.dumps(result.to_dict(), indent=2))


if __name__ == "__main__":
    main()
