from __future__ import annotations

import argparse
import datetime as dt
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

if __package__:
    from .geometry import GeoPoint, coerce_longitude_for_grid, destination_point, great_circle_waypoints, initial_bearing_deg, wrap_longitude
    from .grid import IonosphereGrid, build_pyiri_grid
    from .tracer import PointToPointRayTracer, PyLapImportError, RayTrace
else:
    from geometry import GeoPoint, coerce_longitude_for_grid, destination_point, great_circle_waypoints, initial_bearing_deg, wrap_longitude
    from grid import IonosphereGrid, build_pyiri_grid
    from tracer import PointToPointRayTracer, PyLapImportError, RayTrace


@dataclass(frozen=True)
class GroundToGroundScenario:
    when: dt.datetime
    tx: GeoPoint
    distance_km: float = 500.0
    bearing_deg: float = 90.0
    rx_alt_km: float | None = None
    frequency_mhz: float = 5.0
    f107: float = 120.0
    ap_daily: float = 8.0
    nhops: int = 1
    homing_tolerance_m: float = 500.0


@dataclass(frozen=True)
class GroundToGroundDemoResult:
    scenario: GroundToGroundScenario
    rx: GeoPoint
    ordinary_ray: RayTrace
    extraordinary_ray: RayTrace
    great_circle_distance_km: float
    great_circle_bearing_deg: float
    ordinary_max_path_alt_km: float
    extraordinary_max_path_alt_km: float

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
                "ground_distance_km": self.great_circle_distance_km,
                "initial_bearing_deg": self.great_circle_bearing_deg,
                "frequency_mhz": self.scenario.frequency_mhz,
                "nhops": self.scenario.nhops,
            },
            "solutions": {
                "o_mode": {
                    **self.ordinary_ray.to_dict(),
                    "ox_mode": 1,
                    "max_path_alt_km": self.ordinary_max_path_alt_km,
                },
                "x_mode": {
                    **self.extraordinary_ray.to_dict(),
                    "ox_mode": -1,
                    "max_path_alt_km": self.extraordinary_max_path_alt_km,
                },
            },
        }


@dataclass(frozen=True)
class ElectronDensitySwath:
    along_track_km: np.ndarray
    altitudes_km: np.ndarray
    electron_density_m3: np.ndarray
    ordinary_ray_along_track_km: np.ndarray
    ordinary_ray_altitudes_km: np.ndarray
    extraordinary_ray_along_track_km: np.ndarray
    extraordinary_ray_altitudes_km: np.ndarray


def default_scenario() -> GroundToGroundScenario:
    return GroundToGroundScenario(
        when=dt.datetime(2020, 1, 15, 0, 0, 0),
        tx=GeoPoint(34.873, -106.614, 1.6),
    )


def _parse_time(text: str) -> dt.datetime:
    when = dt.datetime.fromisoformat(text)
    if when.tzinfo is not None:
        when = when.astimezone(dt.timezone.utc).replace(tzinfo=None)
    return when


def _surface_distance_km(start: GeoPoint, end: GeoPoint) -> float:
    lat1 = math.radians(start.lat_deg)
    lat2 = math.radians(end.lat_deg)
    dlat = lat2 - lat1
    dlon = math.radians(((end.lon_deg - start.lon_deg + 180.0) % 360.0) - 180.0)
    a = math.sin(dlat / 2.0) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2
    return 6371.0088 * 2.0 * math.asin(min(1.0, math.sqrt(max(a, 0.0))))


def _raytrace_kwargs(scenario: GroundToGroundScenario) -> dict[str, float | int | str]:
    return {
        "f107": scenario.f107,
        "ap_daily": scenario.ap_daily,
        "nhops": scenario.nhops,
        "homing_tolerance_m": scenario.homing_tolerance_m,
        "alt_min_km": 60.0,
        "alt_max_km": 500.0,
        "alt_step_km": 5.0,
        "lat_step_deg": 0.5,
        "lon_step_deg": 0.5,
        "lat_margin_deg": 3.0,
        "lon_margin_deg": 3.0,
        "d_region_model": "fpt2018",
    }


def build_demo_grid(scenario: GroundToGroundScenario, rx: GeoPoint) -> IonosphereGrid:
    kwargs = _raytrace_kwargs(scenario)
    return build_pyiri_grid(
        scenario.when,
        scenario.tx,
        rx,
        f107=float(kwargs["f107"]),
        ap_daily=float(kwargs["ap_daily"]),
        alt_min_km=float(kwargs["alt_min_km"]),
        alt_max_km=float(kwargs["alt_max_km"]),
        alt_step_km=float(kwargs["alt_step_km"]),
        lat_step_deg=float(kwargs["lat_step_deg"]),
        lon_step_deg=float(kwargs["lon_step_deg"]),
        lat_margin_deg=float(kwargs["lat_margin_deg"]),
        lon_margin_deg=float(kwargs["lon_margin_deg"]),
        d_region_model=str(kwargs["d_region_model"]),
    )


def _ray_along_track_km(tx: GeoPoint, ray_path: dict) -> tuple[np.ndarray, np.ndarray]:
    ray_lats = np.asarray(ray_path.get("lat", []), dtype=float)
    ray_lons = np.asarray(ray_path.get("lon", []), dtype=float)
    ray_alts = np.asarray(ray_path.get("height", []), dtype=float)
    valid = np.isfinite(ray_lats) & np.isfinite(ray_lons) & np.isfinite(ray_alts) & (ray_alts < 1e40)
    ray_lats = ray_lats[valid]
    ray_lons = np.asarray([wrap_longitude(float(value)) for value in ray_lons[valid]], dtype=float)
    ray_alts = ray_alts[valid]
    along_track_km = np.asarray(
        [
            _surface_distance_km(tx, GeoPoint(float(lat), float(lon), 0.0))
            for lat, lon in zip(ray_lats, ray_lons)
        ],
        dtype=float,
    )
    return np.maximum.accumulate(along_track_km), ray_alts


def _max_path_alt_km(ray: RayTrace) -> float:
    heights = np.asarray(ray.path.get("height", []), dtype=float)
    heights = heights[np.isfinite(heights) & (heights < 1e40)]
    return float(np.max(heights)) if heights.size else float("nan")


def build_electron_density_swath(
    result: GroundToGroundDemoResult,
    *,
    waypoint_count: int = 256,
) -> ElectronDensitySwath:
    tx_ground = GeoPoint(result.scenario.tx.lat_deg, result.scenario.tx.lon_deg, 0.0)
    rx_ground = GeoPoint(result.rx.lat_deg, result.rx.lon_deg, 0.0)
    grid = build_demo_grid(result.scenario, result.rx)
    track_lats_deg, track_lons_deg = great_circle_waypoints(tx_ground, rx_ground, count=waypoint_count)
    along_track_km = np.linspace(0.0, result.great_circle_distance_km, waypoint_count, dtype=float)
    sample_lons_deg = np.asarray(
        [coerce_longitude_for_grid(float(lon), grid.longitudes_deg) for lon in track_lons_deg],
        dtype=float,
    )
    lat_mesh, alt_mesh = np.meshgrid(np.asarray(track_lats_deg, dtype=float), grid.altitudes_km, indexing="xy")
    lon_mesh, _ = np.meshgrid(sample_lons_deg, grid.altitudes_km, indexing="xy")
    points = np.column_stack((lat_mesh.ravel(), lon_mesh.ravel(), alt_mesh.ravel()))
    density_interp = RegularGridInterpolator(
        (grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km),
        np.asarray(grid.iono_en_grid, dtype=float) * 1e6,
        bounds_error=False,
        fill_value=np.nan,
    )
    density_m3 = density_interp(points).reshape(grid.altitudes_km.size, along_track_km.size)
    ordinary_ray_along_track_km, ordinary_ray_altitudes_km = _ray_along_track_km(tx_ground, result.ordinary_ray.path)
    extraordinary_ray_along_track_km, extraordinary_ray_altitudes_km = _ray_along_track_km(tx_ground, result.extraordinary_ray.path)
    return ElectronDensitySwath(
        along_track_km=along_track_km,
        altitudes_km=np.asarray(grid.altitudes_km, dtype=float),
        electron_density_m3=density_m3,
        ordinary_ray_along_track_km=ordinary_ray_along_track_km,
        ordinary_ray_altitudes_km=ordinary_ray_altitudes_km,
        extraordinary_ray_along_track_km=extraordinary_ray_along_track_km,
        extraordinary_ray_altitudes_km=extraordinary_ray_altitudes_km,
    )


def _centers_to_edges(values: np.ndarray) -> np.ndarray:
    centers = np.asarray(values, dtype=float)
    if centers.size == 1:
        return np.array([centers[0] - 0.5, centers[0] + 0.5], dtype=float)
    edges = np.empty(centers.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def plot_altitude_profile(path: Path, result: GroundToGroundDemoResult, swath: ElectronDensitySwath | None = None) -> Path:
    import matplotlib.pyplot as plt  # type: ignore

    swath = swath or build_electron_density_swath(result)
    along_edges = _centers_to_edges(swath.along_track_km)
    alt_edges = _centers_to_edges(swath.altitudes_km)

    # Treat undefined cells as zero density so the lower/background region is
    # still visible instead of dropping out of the pcolormesh.
    density_display_m3 = np.where(
        np.isfinite(swath.electron_density_m3) & (swath.electron_density_m3 > 0.0),
        swath.electron_density_m3,
        0.0,
    )
    plasma_frequency_mhz = 8.98e-6 * np.sqrt(np.maximum(density_display_m3, 0.0))
    finite_positive = plasma_frequency_mhz[plasma_frequency_mhz > 0.0]
    vmin = 0.0
    vmax = float(np.max(finite_positive)) if finite_positive.size else 10.0

    # The modeled ionosphere starts at 60 km, but for display we want the full
    # region below that altitude to appear as zero-density background instead of
    # the axes face color.
    plotted_plasma_frequency_mhz = plasma_frequency_mhz
    plotted_alt_edges = alt_edges
    if alt_edges.size >= 2 and alt_edges[0] > 0.0:
        zero_row = np.zeros((1, plasma_frequency_mhz.shape[1]), dtype=float)
        plotted_plasma_frequency_mhz = np.vstack((zero_row, plasma_frequency_mhz))
        plotted_alt_edges = np.concatenate(([0.0], alt_edges))

    ray_top_km = max(
        float(np.nanmax(swath.ordinary_ray_altitudes_km)),
        float(np.nanmax(swath.extraordinary_ray_altitudes_km)),
    )
    display_xmin_km = 0.0
    display_xmax_km = result.great_circle_distance_km
    display_ymin_km = 0.0
    display_ymax_km = max(220.0, ray_top_km + 20.0)

    row_visible = (plotted_alt_edges[:-1] < display_ymax_km) & (plotted_alt_edges[1:] > display_ymin_km)
    col_visible = (along_edges[:-1] < display_xmax_km) & (along_edges[1:] > display_xmin_km)
    visible_values = plotted_plasma_frequency_mhz[np.ix_(row_visible, col_visible)]
    finite_visible = visible_values[np.isfinite(visible_values)]
    if finite_visible.size:
        vmax = float(np.max(finite_visible))
        if vmax <= vmin:
            vmax = vmin + 1e-6

    along_mesh, alt_mesh = np.meshgrid(along_edges, plotted_alt_edges, indexing="xy")

    fig, ax = plt.subplots(figsize=(9.6, 6.2))
    surface = ax.pcolormesh(
        along_mesh,
        alt_mesh,
        plotted_plasma_frequency_mhz,
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
        shading="flat",
        antialiased=False,
        rasterized=True,
    )
    ax.plot(
        swath.ordinary_ray_along_track_km,
        swath.ordinary_ray_altitudes_km,
        color="white",
        linewidth=2.8,
        label="O-mode ray",
    )
    ax.plot(
        swath.extraordinary_ray_along_track_km,
        swath.extraordinary_ray_altitudes_km,
        color="tab:orange",
        linewidth=2.4,
        label="X-mode ray",
    )
    ax.scatter(
        [0.0, result.great_circle_distance_km],
        [result.scenario.tx.alt_km, result.rx.alt_km],
        c=["gold", "tab:red"],
        edgecolors="black",
        linewidths=0.8,
        s=70,
        zorder=5,
    )
    ordinary_apex_idx = int(np.argmax(swath.ordinary_ray_altitudes_km))
    extraordinary_apex_idx = int(np.argmax(swath.extraordinary_ray_altitudes_km))
    ax.scatter(
        [swath.ordinary_ray_along_track_km[ordinary_apex_idx]],
        [swath.ordinary_ray_altitudes_km[ordinary_apex_idx]],
        color="white",
        s=35,
        zorder=6,
    )
    ax.scatter(
        [swath.extraordinary_ray_along_track_km[extraordinary_apex_idx]],
        [swath.extraordinary_ray_altitudes_km[extraordinary_apex_idx]],
        color="tab:orange",
        edgecolors="black",
        linewidths=0.6,
        s=35,
        zorder=6,
    )
    ax.set_xlim(0.0, result.great_circle_distance_km)
    ax.set_ylim(display_ymin_km, display_ymax_km)
    ax.set_xlabel("Ground-track distance from transmitter (km)")
    ax.set_ylabel("Altitude (km)")
    ax.set_title(
        f"Ground-to-Ground {result.scenario.frequency_mhz:.1f} MHz Raytrace: Electron Density Swath with O/X Rays"
    )
    ax.set_axisbelow(True)
    ax.minorticks_on()
    ax.grid(True, which="major", alpha=0.28, linewidth=0.7)
    ax.grid(True, which="minor", alpha=0.12, linewidth=0.5, linestyle=":")
    ax.legend(loc="upper right", framealpha=0.95)
    subtitle = (
        f"O: elev={result.ordinary_ray.launch_elevation_deg:.2f} deg, miss={result.ordinary_ray.error_m:.1f} m, "
        f"abs={result.ordinary_ray.total_absorption_db:.2f} dB | "
        f"X: elev={result.extraordinary_ray.launch_elevation_deg:.2f} deg, miss={result.extraordinary_ray.error_m:.1f} m, "
        f"abs={result.extraordinary_ray.total_absorption_db:.2f} dB"
    )
    ax.text(
        0.01,
        0.02,
        subtitle,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.35, "pad": 4, "edgecolor": "none"},
    )
    cbar = fig.colorbar(surface, ax=ax, pad=0.02, label="Plasma frequency (MHz), undefined=0")
    cbar.ax.tick_params(labelsize=10)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    return path


def run_ground_to_ground_demo(scenario: GroundToGroundScenario | None = None) -> GroundToGroundDemoResult:
    scenario = scenario or default_scenario()
    rx = destination_point(
        scenario.tx,
        scenario.bearing_deg,
        scenario.distance_km,
        alt_km=scenario.tx.alt_km if scenario.rx_alt_km is None else scenario.rx_alt_km,
    )
    tracer = PointToPointRayTracer()
    grid = build_demo_grid(scenario, rx)
    ordinary_ray = tracer.trace_link(
        tx=scenario.tx,
        rx=rx,
        frequency_mhz=scenario.frequency_mhz,
        grid=grid,
        ox_mode=1,
        nhops=scenario.nhops,
        homing_tolerance_m=scenario.homing_tolerance_m,
    )
    extraordinary_ray = tracer.trace_link(
        tx=scenario.tx,
        rx=rx,
        frequency_mhz=scenario.frequency_mhz,
        grid=grid,
        ox_mode=-1,
        nhops=scenario.nhops,
        homing_tolerance_m=scenario.homing_tolerance_m,
    )
    return GroundToGroundDemoResult(
        scenario=scenario,
        rx=rx,
        ordinary_ray=ordinary_ray,
        extraordinary_ray=extraordinary_ray,
        great_circle_distance_km=_surface_distance_km(scenario.tx, rx),
        great_circle_bearing_deg=initial_bearing_deg(scenario.tx, rx),
        ordinary_max_path_alt_km=_max_path_alt_km(ordinary_ray),
        extraordinary_max_path_alt_km=_max_path_alt_km(extraordinary_ray),
    )


def build_parser() -> argparse.ArgumentParser:
    scenario = default_scenario()
    parser = argparse.ArgumentParser(
        description="Ground-to-ground 3D homing example using the real PyIRI + pylap raytracer.",
    )
    parser.add_argument("--time", default=scenario.when.isoformat(), help="UTC time in ISO-8601 format.")
    parser.add_argument("--tx", nargs=3, type=float, metavar=("LAT", "LON", "ALT_KM"),
                        default=(scenario.tx.lat_deg, scenario.tx.lon_deg, scenario.tx.alt_km))
    parser.add_argument("--distance-km", type=float, default=scenario.distance_km)
    parser.add_argument("--bearing-deg", type=float, default=scenario.bearing_deg)
    parser.add_argument("--rx-alt-km", type=float, default=scenario.tx.alt_km)
    parser.add_argument("--freq-mhz", type=float, default=scenario.frequency_mhz)
    parser.add_argument("--f107", type=float, default=scenario.f107)
    parser.add_argument("--ap-daily", type=float, default=scenario.ap_daily)
    parser.add_argument("--nhops", type=int, default=scenario.nhops)
    parser.add_argument("--homing-tolerance-m", type=float, default=scenario.homing_tolerance_m)
    parser.add_argument("--profile-out", type=Path, default=None, help="Optional PNG output path for the electron-density swath altitude profile.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    scenario = GroundToGroundScenario(
        when=_parse_time(args.time),
        tx=GeoPoint(*args.tx),
        distance_km=float(args.distance_km),
        bearing_deg=float(args.bearing_deg),
        rx_alt_km=float(args.rx_alt_km),
        frequency_mhz=float(args.freq_mhz),
        f107=float(args.f107),
        ap_daily=float(args.ap_daily),
        nhops=int(args.nhops),
        homing_tolerance_m=float(args.homing_tolerance_m),
    )
    try:
        result = run_ground_to_ground_demo(scenario)
    except PyLapImportError as exc:
        raise SystemExit(str(exc)) from exc
    if args.profile_out is not None:
        plot_altitude_profile(args.profile_out.expanduser(), result)
    print(json.dumps(result.to_dict(), indent=2))


if __name__ == "__main__":
    main()
