from __future__ import annotations

import argparse
import csv
import datetime as dt
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Sequence

import numpy as np
from sgp4.api import Satrec

if __package__:
    from .absorption import effective_collision_frequency
    from .geometry import GeoPoint, WGS84_A_M, WGS84_E2, great_circle_waypoints, initial_bearing_deg, llh_to_ecef, ray_point_distance, wrap_longitude
    from .grid import IonosphereGrid, build_pyiri_grid, build_pyiri_grid_from_axes
    from .indices import resolve_space_weather_indices
    from .tracer import PointToPointRayTracer, PyLapImportError, RayTrace
else:
    from absorption import effective_collision_frequency
    from geometry import GeoPoint, WGS84_A_M, WGS84_E2, great_circle_waypoints, initial_bearing_deg, llh_to_ecef, ray_point_distance, wrap_longitude
    from grid import IonosphereGrid, build_pyiri_grid, build_pyiri_grid_from_axes
    from indices import resolve_space_weather_indices
    from tracer import PointToPointRayTracer, PyLapImportError, RayTrace


C_KM_PER_S = 299_792.458
C_M_PER_S = 299_792_458.0
DEFAULT_TLE_FILE = Path("~/superdarn/digital_rf_tools/artifacts/iss_stations.tle").expanduser()
DEFAULT_MEASUREMENT_FILE = Path("~/superdarn/digital_rf_tools/tmp/pulsed_meteor_radar_detection_METnwDEU_Collm_Skiymet_36.2000MHz.csv").expanduser()
DEFAULT_START_UTC = dt.datetime(2026, 6, 1, 15, 58, 35, tzinfo=dt.timezone.utc)


@dataclass(frozen=True)
class CollmSkiymetRadar:
    network: str = "METnwDEU"
    site: str = "Collm"
    system: str = "Skiymet"
    tx: GeoPoint = GeoPoint(51.309, 13.003, 0.170)
    frequency_mhz: float = 36.2000
    prf_hz: float = 625.0
    chip_rate_hz: float = 100_000.0
    code: str = "Barker-7"
    channel_rate_hz: float = 400_000.0

    @property
    def freq_hz(self) -> float:
        return self.frequency_mhz * 1e6

    @property
    def unambiguous_range_km(self) -> float:
        return C_KM_PER_S / self.prf_hz

    @property
    def doppler_nyquist_hz(self) -> float:
        return 0.5 * self.prf_hz

    @property
    def chip_samples(self) -> int:
        return max(1, int(round(self.channel_rate_hz / self.chip_rate_hz)))

    @property
    def pulse_samples(self) -> int:
        return 7 * self.chip_samples

    @property
    def label(self) -> str:
        return f"{self.network}/{self.site} ({self.system})"

    @property
    def title(self) -> str:
        return (
            f"{self.network}/{self.site} ({self.system}): "
            f"{self.code} ({self.chip_samples} samples/chip, {self.pulse_samples} samples/pulse)"
        )


@dataclass(frozen=True)
class SNRModel:
    reference_snr_db: float = 8.8
    reference_range_km: float = 700.0
    reference_elevation_deg: float = 35.0
    range_power_exponent: float = 5.0
    beam_power_exponent: float = 3.0
    ripple_amplitude_db: float = 0.9
    ripple_period_km: float = 115.0
    ripple_phase_rad: float = 0.6
    floor_snr_db: float = -2.0


@dataclass(frozen=True)
class HookeWaveModel:
    amplitude_fraction: float = 0.16
    horizontal_wavelength_km: float = 82.0
    bearing_deg: float | None = None
    bearing_offset_deg: float = 0.0
    period_seconds: float = 12.6
    phase_rad: float = 0.75
    vertical_center_km: float = 285.0
    vertical_sigma_km: float = 70.0
    snr_coupling_db: float = 1.45


@dataclass(frozen=True)
class RayTubeModel:
    delta_elevation_deg: float = 0.2
    delta_bearing_deg: float = 0.2


@dataclass(frozen=True)
class TimeSeriesConfig:
    tle_file: Path = DEFAULT_TLE_FILE
    start_utc: dt.datetime = DEFAULT_START_UTC
    seconds: float = 90.0
    step_seconds: float = 1.0
    ephemeris_time_shift_seconds: float = 0.0
    min_elevation_deg: float = 0.0
    range_offset_km: float = 0.0
    doppler_offset_hz: float = 0.0
    homing_tolerance_m: float = 500.0
    f107: float | None = None
    f107a: float | None = None
    snr_model: SNRModel = SNRModel()
    hooke_wave: HookeWaveModel = HookeWaveModel()
    ray_tube_model: RayTubeModel = RayTubeModel()
    real_grid_lat_step_deg: float = 0.5
    real_grid_lon_step_deg: float = 0.5
    real_grid_lat_margin_deg: float = 3.0
    real_grid_lon_margin_deg: float = 3.0
    real_grid_alt_min_km: float = 60.0
    real_grid_alt_max_km: float = 700.0
    real_grid_alt_step_km: float = 5.0
    real_d_region_model: str = "fpt2018"
    noise_seed: int = 0


@dataclass(frozen=True)
class PredictionPoint:
    time_utc: dt.datetime
    satellite: GeoPoint
    los_bearing_deg: float
    los_elevation_deg: float
    absolute_range_km: float
    folded_range_km: float
    absolute_doppler_hz: float
    aliased_doppler_hz: float
    predicted_peak_snr_db: float
    launch_bearing_deg: float | None
    launch_elevation_deg: float | None
    group_range_km: float | None
    total_absorption_db: float | None
    ray_tube_area_m2: float | None
    ray_tube_gain_db: float | None
    max_path_alt_km: float | None
    home: bool
    error_m: float | None


@dataclass(frozen=True)
class FoF2MapSnapshot:
    time_utc: dt.datetime
    latitudes_deg: np.ndarray
    longitudes_deg: np.ndarray
    fof2_mhz: np.ndarray
    background_fof2_mhz: np.ndarray
    wave_origin_lat_deg: float
    wave_origin_lon_deg: float
    wave_bearing_deg: float
    wave_elapsed_seconds: float
    hooke_wave: HookeWaveModel
    track_latitudes_deg: np.ndarray
    track_longitudes_deg: np.ndarray
    transmitter_lat_deg: float
    transmitter_lon_deg: float
    cpa_lat_deg: float
    cpa_lon_deg: float
    cpa_index: int
    f107: float
    f107a: float
    display_lon_min_deg: float
    display_lon_max_deg: float
    display_lat_min_deg: float
    display_lat_max_deg: float


@dataclass(frozen=True)
class PredictionSeries:
    radar: CollmSkiymetRadar
    config: TimeSeriesConfig
    points: tuple[PredictionPoint, ...]
    range_axis_km: np.ndarray
    image_snr_db: np.ndarray


@dataclass(frozen=True)
class MeasurementPoint:
    time_utc: dt.datetime
    peak_range_km: float
    peak_snr_db: float
    doppler_hz: float


@dataclass(frozen=True)
class SwathSnapshot:
    time_utc: dt.datetime
    along_track_km: np.ndarray
    altitudes_km: np.ndarray
    plasma_freq_mhz: np.ndarray
    ray_along_track_km: np.ndarray
    ray_altitudes_km: np.ndarray
    tx_ground_label: str
    rx_ground_label: str
    surface_range_km: float
    target_alt_km: float
    f107: float
    f107a: float


@dataclass(frozen=True)
class RayTubeMetrics:
    area_m2_per_rad2: float
    gain_db: float


@dataclass(frozen=True)
class RealPassContext:
    tracer: PointToPointRayTracer
    background_grid: IonosphereGrid
    phase_coordinate_km: np.ndarray


@dataclass(frozen=True)
class GeometryFitResult:
    ephemeris_time_shift_seconds: float
    range_offset_km: float
    doppler_offset_hz: float
    rms_range_km: float
    rms_doppler_hz: float


@dataclass(frozen=True)
class SNRFitResult:
    offset_db: float
    scale: float
    rms_db: float


def _parse_time(text: str) -> dt.datetime:
    when = dt.datetime.fromisoformat(text)
    if when.tzinfo is None:
        return when.replace(tzinfo=dt.timezone.utc)
    return when.astimezone(dt.timezone.utc)


def load_tle(path: Path) -> Satrec:
    lines = [line.rstrip() for line in path.read_text().splitlines() if line.strip()]
    if len(lines) < 2:
        raise RuntimeError(f"TLE file {path} does not contain two lines.")
    if lines[0].startswith("1 "):
        line1, line2 = lines[0], lines[1]
    else:
        line1, line2 = lines[1], lines[2]
    return Satrec.twoline2rv(line1, line2)


def load_measurements(path: Path, *, start_utc: dt.datetime | None = None, end_utc: dt.datetime | None = None) -> tuple[MeasurementPoint, ...]:
    points: list[MeasurementPoint] = []
    with path.expanduser().open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            when = _parse_time(row["utc_time"])
            if start_utc is not None and when < start_utc:
                continue
            if end_utc is not None and when > end_utc:
                continue
            peak_range_km_text = row.get("peak_range_km", "")
            peak_snr_db_text = row.get("peak_snr_db", "")
            doppler_text = row.get("predicted_doppler_hz", "") or row.get("residual_doppler_hz", "")
            if not peak_range_km_text or not peak_snr_db_text or not doppler_text:
                continue
            peak_range_km = float(peak_range_km_text)
            peak_snr_db = float(peak_snr_db_text)
            doppler_hz = float(doppler_text)
            if not (math.isfinite(peak_range_km) and math.isfinite(peak_snr_db) and math.isfinite(doppler_hz)):
                continue
            points.append(
                MeasurementPoint(
                    time_utc=when,
                    peak_range_km=peak_range_km,
                    peak_snr_db=peak_snr_db,
                    doppler_hz=doppler_hz,
                )
            )
    return tuple(points)


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


def _ecef_to_teme(r_m: np.ndarray, jd_ut1: float) -> np.ndarray:
    theta = _gmst_from_jd(jd_ut1)
    c = math.cos(theta)
    s = math.sin(theta)
    r_m = np.asarray(r_m, dtype=np.float64)
    return np.array(
        [
            c * r_m[0] - s * r_m[1],
            s * r_m[0] + c * r_m[1],
            r_m[2],
        ],
        dtype=np.float64,
    )


def _sat_teme_at_unix(sat: Satrec, unix_time_s: float) -> np.ndarray:
    jd = unix_time_s / 86400.0 + 2440587.5
    jd0 = math.floor(jd)
    fr = jd - jd0
    err, r_km, _v_km_s = sat.sgp4(jd0, fr)
    if err != 0:
        raise RuntimeError(f"SGP4 propagation failed with code {err} at unix time {unix_time_s}.")
    return np.asarray(r_km, dtype=np.float64)


def _sat_ecef_at_unix(sat: Satrec, unix_time_s: float) -> np.ndarray:
    jd = unix_time_s / 86400.0 + 2440587.5
    return _teme_to_ecef(_sat_teme_at_unix(sat, unix_time_s), jd)


def _ecef_to_geo(point_xyz_m: np.ndarray) -> GeoPoint:
    point = np.asarray(point_xyz_m, dtype=np.float64)
    x = float(point[0])
    y = float(point[1])
    z = float(point[2])
    lon = math.atan2(y, x)
    p = math.hypot(x, y)
    lat = math.atan2(z, p * (1.0 - WGS84_E2))
    alt_m = 0.0
    for _ in range(7):
        sin_lat = math.sin(lat)
        cos_lat = math.cos(lat)
        prime_vertical = WGS84_A_M / math.sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat)
        safe_cos = cos_lat if abs(cos_lat) > 1e-12 else math.copysign(1e-12, cos_lat if cos_lat != 0.0 else 1.0)
        alt_m = p / safe_cos - prime_vertical
        denom = max(prime_vertical + alt_m, 1.0)
        lat = math.atan2(z, p * (1.0 - WGS84_E2 * prime_vertical / denom))
    return GeoPoint(math.degrees(lat), math.degrees(lon), alt_m / 1000.0)


def _one_way_light_time_s(sat: Satrec, receive_unix_time_s: float, radar_ecef_m: np.ndarray, max_iter: int = 4) -> float:
    jd_rx = receive_unix_time_s / 86400.0 + 2440587.5
    radar_teme_m = _ecef_to_teme(radar_ecef_m, jd_rx)
    tau_s = 0.0
    for _ in range(max(int(max_iter), 1)):
        transmit_unix_time_s = receive_unix_time_s - tau_s
        sat_teme_m = _sat_teme_at_unix(sat, transmit_unix_time_s) * 1000.0
        rho_m = float(np.linalg.norm(sat_teme_m - radar_teme_m))
        new_tau_s = rho_m / C_M_PER_S
        if abs(new_tau_s - tau_s) < 1e-12:
            tau_s = new_tau_s
            break
        tau_s = new_tau_s
    return tau_s


def predict_delay_doppler(
    sat: Satrec,
    frame_times_s: np.ndarray,
    radar_ecef_m: np.ndarray,
    freq_hz: float,
) -> tuple[np.ndarray, np.ndarray]:
    delays_s = np.zeros(frame_times_s.size, dtype=np.float64)
    dopplers_hz = np.zeros(frame_times_s.size, dtype=np.float64)
    dt_seconds = 0.5
    for idx, ts in enumerate(frame_times_s):
        r0_m = _one_way_light_time_s(sat, float(ts), radar_ecef_m) * C_M_PER_S
        rm_m = _one_way_light_time_s(sat, float(ts) - dt_seconds, radar_ecef_m) * C_M_PER_S
        rp_m = _one_way_light_time_s(sat, float(ts) + dt_seconds, radar_ecef_m) * C_M_PER_S
        range_rate_mps = (rp_m - rm_m) / (2.0 * dt_seconds)
        delays_s[idx] = r0_m / C_M_PER_S
        dopplers_hz[idx] = -freq_hz * range_rate_mps / C_M_PER_S
    return delays_s, dopplers_hz


def _line_of_sight_angles(tx: GeoPoint, rx: GeoPoint) -> tuple[float, float, float]:
    tx_xyz = np.asarray(llh_to_ecef(tx.lat_deg, tx.lon_deg, tx.alt_km * 1000.0), dtype=np.float64)
    rx_xyz = np.asarray(llh_to_ecef(rx.lat_deg, rx.lon_deg, rx.alt_km * 1000.0), dtype=np.float64)
    if tx_xyz.ndim > 1:
        tx_xyz = tx_xyz[0]
    if rx_xyz.ndim > 1:
        rx_xyz = rx_xyz[0]
    los = rx_xyz - tx_xyz
    slant_range_km = float(np.linalg.norm(los) / 1000.0)

    lat_rad = math.radians(tx.lat_deg)
    lon_rad = math.radians(tx.lon_deg)
    east = np.array([-math.sin(lon_rad), math.cos(lon_rad), 0.0], dtype=np.float64)
    north = np.array(
        [
            -math.sin(lat_rad) * math.cos(lon_rad),
            -math.sin(lat_rad) * math.sin(lon_rad),
            math.cos(lat_rad),
        ],
        dtype=np.float64,
    )
    up = np.array(
        [
            math.cos(lat_rad) * math.cos(lon_rad),
            math.cos(lat_rad) * math.sin(lon_rad),
            math.sin(lat_rad),
        ],
        dtype=np.float64,
    )
    east_component = float(np.dot(los, east))
    north_component = float(np.dot(los, north))
    up_component = float(np.dot(los, up))

    bearing_deg = math.degrees(math.atan2(east_component, north_component))
    elevation_deg = math.degrees(math.atan2(up_component, math.hypot(east_component, north_component)))
    return bearing_deg, elevation_deg, slant_range_km


def _alias_to_prf_band(values_hz: np.ndarray, prf_hz: float) -> np.ndarray:
    nyquist = 0.5 * prf_hz
    return ((values_hz + nyquist) % prf_hz) - nyquist


def _wrapped_range_residual_km(predicted_km: np.ndarray, measured_km: np.ndarray, unambiguous_range_km: float) -> np.ndarray:
    residual = np.asarray(predicted_km, dtype=np.float64) - np.asarray(measured_km, dtype=np.float64)
    half = 0.5 * float(unambiguous_range_km)
    return ((residual + half) % float(unambiguous_range_km)) - half


def _fit_geometry_to_measurements(
    radar: CollmSkiymetRadar,
    tle_file: Path,
    measurements: Sequence[MeasurementPoint],
    *,
    initial_ephemeris_time_shift_seconds: float = 0.0,
    initial_range_offset_km: float = 0.0,
    initial_doppler_offset_hz: float = 0.0,
) -> GeometryFitResult:
    from scipy.optimize import minimize

    if not measurements:
        raise ValueError("at least one measurement is required for geometry fitting")

    sat = load_tle(tle_file.expanduser())
    times = np.asarray([point.time_utc.timestamp() for point in measurements], dtype=np.float64)
    measured_ranges_km = np.asarray([point.peak_range_km for point in measurements], dtype=np.float64)
    measured_dopplers_hz = np.asarray([point.doppler_hz for point in measurements], dtype=np.float64)
    radar_ecef_m = np.asarray(llh_to_ecef(radar.tx.lat_deg, radar.tx.lon_deg, radar.tx.alt_km * 1000.0), dtype=np.float64)
    if radar_ecef_m.ndim > 1:
        radar_ecef_m = radar_ecef_m[0]

    def objective(params: np.ndarray) -> float:
        shifted_times = times + float(params[0])
        delays_s, absolute_dopplers_hz = predict_delay_doppler(sat, shifted_times, radar_ecef_m, radar.freq_hz)
        folded_ranges_km = np.mod(delays_s * C_KM_PER_S + float(params[1]), radar.unambiguous_range_km)
        aliased_dopplers_hz = _alias_to_prf_band(absolute_dopplers_hz + float(params[2]), radar.prf_hz)
        range_residual_km = _wrapped_range_residual_km(folded_ranges_km, measured_ranges_km, radar.unambiguous_range_km)
        doppler_residual_hz = aliased_dopplers_hz - measured_dopplers_hz
        return float(np.mean((range_residual_km / 1.5) ** 2 + (doppler_residual_hz / 15.0) ** 2))

    x0 = np.array(
        [
            initial_ephemeris_time_shift_seconds,
            initial_range_offset_km,
            initial_doppler_offset_hz,
        ],
        dtype=np.float64,
    )
    result = minimize(
        objective,
        x0,
        method="Nelder-Mead",
        options={"xatol": 1e-3, "fatol": 1e-4, "maxfev": 300, "disp": False},
    )
    best = np.asarray(result.x, dtype=np.float64)
    shifted_times = times + float(best[0])
    delays_s, absolute_dopplers_hz = predict_delay_doppler(sat, shifted_times, radar_ecef_m, radar.freq_hz)
    folded_ranges_km = np.mod(delays_s * C_KM_PER_S + float(best[1]), radar.unambiguous_range_km)
    aliased_dopplers_hz = _alias_to_prf_band(absolute_dopplers_hz + float(best[2]), radar.prf_hz)
    range_residual_km = _wrapped_range_residual_km(folded_ranges_km, measured_ranges_km, radar.unambiguous_range_km)
    doppler_residual_hz = aliased_dopplers_hz - measured_dopplers_hz
    return GeometryFitResult(
        ephemeris_time_shift_seconds=float(best[0]),
        range_offset_km=float(best[1]),
        doppler_offset_hz=float(best[2]),
        rms_range_km=float(np.sqrt(np.mean(range_residual_km ** 2))),
        rms_doppler_hz=float(np.sqrt(np.mean(doppler_residual_hz ** 2))),
    )


def _fit_snr_to_measurements(series: PredictionSeries, measurements: Sequence[MeasurementPoint]) -> SNRFitResult:
    measurement_by_time = {point.time_utc: point for point in measurements}
    modeled: list[float] = []
    observed: list[float] = []
    for point in series.points:
        measurement = measurement_by_time.get(point.time_utc)
        if measurement is None or not math.isfinite(point.predicted_peak_snr_db):
            continue
        modeled.append(float(point.predicted_peak_snr_db))
        observed.append(float(measurement.peak_snr_db))
    if not modeled:
        raise ValueError("no overlapping modeled/measured SNR samples for SNR fitting")
    x = np.asarray(modeled, dtype=np.float64)
    y = np.asarray(observed, dtype=np.float64)
    design = np.column_stack((np.ones_like(x), x))
    coeffs, *_ = np.linalg.lstsq(design, y, rcond=None)
    fitted = design @ coeffs
    rms_db = float(np.sqrt(np.mean((fitted - y) ** 2)))
    return SNRFitResult(
        offset_db=float(coeffs[0]),
        scale=float(coeffs[1]),
        rms_db=rms_db,
    )


def _apply_snr_fit(series: PredictionSeries, fit: SNRFitResult) -> PredictionSeries:
    calibrated_points: list[PredictionPoint] = []
    for point in series.points:
        calibrated_snr = point.predicted_peak_snr_db
        if math.isfinite(calibrated_snr):
            calibrated_snr = fit.offset_db + fit.scale * calibrated_snr
        calibrated_points.append(
            PredictionPoint(
                time_utc=point.time_utc,
                satellite=point.satellite,
                los_bearing_deg=point.los_bearing_deg,
                los_elevation_deg=point.los_elevation_deg,
                absolute_range_km=point.absolute_range_km,
                folded_range_km=point.folded_range_km,
                absolute_doppler_hz=point.absolute_doppler_hz,
                aliased_doppler_hz=point.aliased_doppler_hz,
                predicted_peak_snr_db=calibrated_snr,
                launch_bearing_deg=point.launch_bearing_deg,
                launch_elevation_deg=point.launch_elevation_deg,
                group_range_km=point.group_range_km,
                total_absorption_db=point.total_absorption_db,
                ray_tube_area_m2=point.ray_tube_area_m2,
                ray_tube_gain_db=point.ray_tube_gain_db,
                max_path_alt_km=point.max_path_alt_km,
                home=point.home,
                error_m=point.error_m,
            )
        )
    range_axis_km, image_snr_db = _build_synthetic_delay_image(series.radar, calibrated_points, series.config.noise_seed)
    return PredictionSeries(
        radar=series.radar,
        config=series.config,
        points=tuple(calibrated_points),
        range_axis_km=range_axis_km,
        image_snr_db=image_snr_db,
    )


def _time_grid(start_utc: dt.datetime, seconds: float, step_seconds: float) -> list[dt.datetime]:
    count = int(math.floor(max(seconds, 0.0) / max(step_seconds, 1e-9))) + 1
    return [start_utc + dt.timedelta(seconds=idx * step_seconds) for idx in range(count)]


def _centers_to_edges(values: np.ndarray) -> np.ndarray:
    centers = np.asarray(values, dtype=np.float64)
    if centers.size == 0:
        return np.zeros(0, dtype=np.float64)
    if centers.size == 1:
        return np.array([centers[0] - 0.5, centers[0] + 0.5], dtype=np.float64)
    edges = np.empty(centers.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = centers[0] - 0.5 * (centers[1] - centers[0])
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    return edges


def _coerce_longitudes_near_reference(longitudes_deg: np.ndarray, reference_deg: float) -> np.ndarray:
    longitudes = np.asarray(longitudes_deg, dtype=np.float64)
    candidates = np.stack((longitudes - 360.0, longitudes, longitudes + 360.0), axis=1)
    index = np.argmin(np.abs(candidates - float(reference_deg)), axis=1)
    return candidates[np.arange(longitudes.size), index]


def _closest_approach_index(points: Sequence[PredictionPoint]) -> int:
    ranges = np.asarray([point.absolute_range_km for point in points], dtype=np.float64)
    finite = np.isfinite(ranges)
    if not np.any(finite):
        return 0
    finite_indices = np.flatnonzero(finite)
    return int(finite_indices[np.argmin(ranges[finite])])


def _surface_distance_km(start: GeoPoint, end: GeoPoint) -> float:
    lat1 = math.radians(start.lat_deg)
    lat2 = math.radians(end.lat_deg)
    dlat = lat2 - lat1
    dlon = math.radians(wrap_longitude(end.lon_deg - start.lon_deg))
    a = math.sin(dlat / 2.0) ** 2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2.0) ** 2
    return 6371.0088 * 2.0 * math.asin(min(1.0, math.sqrt(max(a, 0.0))))


def _local_east_north_km(origin: GeoPoint, lat_deg: np.ndarray | float, lon_deg: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
    lat = np.asarray(lat_deg, dtype=np.float64)
    lon = np.asarray(lon_deg, dtype=np.float64)
    origin_lat_rad = math.radians(origin.lat_deg)
    dlat_rad = np.deg2rad(lat - origin.lat_deg)
    dlon_rad = np.deg2rad(np.asarray([wrap_longitude(float(value) - origin.lon_deg) for value in np.ravel(lon)], dtype=np.float64)).reshape(lon.shape)
    north_km = 6371.0088 * dlat_rad
    east_km = 6371.0088 * math.cos(origin_lat_rad) * dlon_rad
    return east_km, north_km


def _along_track_coordinate_km(origin: GeoPoint, bearing_deg: float, lat_deg: np.ndarray | float, lon_deg: np.ndarray | float) -> np.ndarray:
    east_km, north_km = _local_east_north_km(origin, lat_deg, lon_deg)
    az = math.radians(bearing_deg)
    return east_km * math.sin(az) + north_km * math.cos(az)


def _resolve_wave_bearing_deg(wave: HookeWaveModel, fallback_bearing_deg: float) -> float:
    base_bearing_deg = fallback_bearing_deg if wave.bearing_deg is None else wave.bearing_deg
    return float(base_bearing_deg + wave.bearing_offset_deg)


def _hooke_wave_fraction(
    wave: HookeWaveModel,
    along_track_km: np.ndarray | float,
    alt_km: np.ndarray | float,
    elapsed_seconds: float,
) -> np.ndarray:
    along = np.asarray(along_track_km, dtype=np.float64)
    alt = np.asarray(alt_km, dtype=np.float64)
    if wave.amplitude_fraction == 0.0:
        return np.zeros(np.broadcast_shapes(along.shape, alt.shape), dtype=np.float64)
    vertical = np.exp(-0.5 * ((alt - wave.vertical_center_km) / max(wave.vertical_sigma_km, 1e-6)) ** 2)
    phase = (
        2.0 * math.pi * (along / max(wave.horizontal_wavelength_km, 1e-6))
        - 2.0 * math.pi * (elapsed_seconds / max(wave.period_seconds, 1e-6))
        + wave.phase_rad
    )
    return wave.amplitude_fraction * vertical * np.cos(phase)


def _apply_hooke_wave_to_density(
    density_m3: np.ndarray,
    along_track_km: np.ndarray,
    altitudes_km: np.ndarray,
    elapsed_seconds: float,
    wave: HookeWaveModel,
) -> np.ndarray:
    density = np.asarray(density_m3, dtype=np.float64)
    along = np.asarray(along_track_km, dtype=np.float64)
    alts = np.asarray(altitudes_km, dtype=np.float64)
    if density.ndim == 2:
        perturbation = _hooke_wave_fraction(wave, along[None, :], alts[:, None], elapsed_seconds)
    elif density.ndim == 3:
        perturbation = _hooke_wave_fraction(wave, along[None, :, :], alts[:, None, None], elapsed_seconds)
    else:
        raise ValueError("density_m3 must be a 2D or 3D array")
    return np.maximum(density * (1.0 + perturbation), 1e-6)


def _build_real_pass_context(
    radar: CollmSkiymetRadar,
    targets: Sequence[GeoPoint],
    background_time_utc: dt.datetime,
    config: TimeSeriesConfig,
    wave_bearing_deg: float,
) -> RealPassContext:
    track_lats = np.asarray([target.lat_deg for target in targets], dtype=np.float64)
    track_lons_wrapped = np.asarray([wrap_longitude(target.lon_deg) for target in targets], dtype=np.float64)
    reference_lon = float(np.mean(track_lons_wrapped)) if track_lons_wrapped.size else float(wrap_longitude(radar.tx.lon_deg))
    track_lons = _coerce_longitudes_near_reference(track_lons_wrapped, reference_lon)

    lat_min = max(-90.0, float(np.min(np.append(track_lats, radar.tx.lat_deg))) - config.real_grid_lat_margin_deg)
    lat_max = min(90.0, float(np.max(np.append(track_lats, radar.tx.lat_deg))) + config.real_grid_lat_margin_deg)
    lon_min = float(np.min(np.append(track_lons, _coerce_longitudes_near_reference(np.array([wrap_longitude(radar.tx.lon_deg)]), reference_lon)[0]))) - config.real_grid_lon_margin_deg
    lon_max = float(np.max(np.append(track_lons, _coerce_longitudes_near_reference(np.array([wrap_longitude(radar.tx.lon_deg)]), reference_lon)[0]))) + config.real_grid_lon_margin_deg

    lat_step = max(config.real_grid_lat_step_deg, 1e-6)
    lon_step = max(config.real_grid_lon_step_deg, 1e-6)
    alt_step = max(config.real_grid_alt_step_km, 1e-6)
    latitudes_deg = np.arange(
        math.floor(lat_min / lat_step) * lat_step,
        math.ceil(lat_max / lat_step) * lat_step + 0.5 * lat_step,
        lat_step,
        dtype=np.float64,
    )
    longitudes_deg = np.arange(
        math.floor(lon_min / lon_step) * lon_step,
        math.ceil(lon_max / lon_step) * lon_step + 0.5 * lon_step,
        lon_step,
        dtype=np.float64,
    )
    alt_max_km = max(config.real_grid_alt_max_km, max(target.alt_km for target in targets) + 50.0)
    altitudes_km = np.arange(
        config.real_grid_alt_min_km,
        alt_max_km + 0.5 * alt_step,
        alt_step,
        dtype=np.float64,
    )
    background_grid = build_pyiri_grid_from_axes(
        background_time_utc.astimezone(dt.timezone.utc).replace(tzinfo=None),
        latitudes_deg,
        longitudes_deg,
        altitudes_km,
        f107=config.f107,
        d_region_model=config.real_d_region_model,
    )
    lon_mesh, lat_mesh = np.meshgrid(background_grid.longitudes_deg, background_grid.latitudes_deg, indexing="xy")
    tx_ground = GeoPoint(radar.tx.lat_deg, radar.tx.lon_deg, 0.0)
    phase_coordinate_km = _along_track_coordinate_km(tx_ground, wave_bearing_deg, lat_mesh, lon_mesh)
    return RealPassContext(
        tracer=PointToPointRayTracer(),
        background_grid=background_grid,
        phase_coordinate_km=phase_coordinate_km,
    )


def _wave_perturbed_grid(
    background_grid: IonosphereGrid,
    phase_coordinate_km: np.ndarray,
    elapsed_seconds: float,
    wave: HookeWaveModel,
) -> IonosphereGrid:
    density_lat_lon_alt_m3 = np.asarray(background_grid.iono_en_grid, dtype=np.float64) * 1e6
    density_alt_lat_lon_m3 = np.transpose(density_lat_lon_alt_m3, (2, 0, 1))
    perturbed_density_alt_lat_lon_m3 = _apply_hooke_wave_to_density(
        density_alt_lat_lon_m3,
        phase_coordinate_km,
        background_grid.altitudes_km,
        elapsed_seconds,
        wave,
    )
    perturbed_density_m3 = np.transpose(perturbed_density_alt_lat_lon_m3, (1, 2, 0))

    collision_freq = np.asarray(background_grid.collision_freq, dtype=np.float64)
    if (
        background_grid.electron_temp_k is not None
        and background_grid.ion_temp_k is not None
        and background_grid.neutral_species_cm3 is not None
    ):
        collision_freq = effective_collision_frequency(
            background_grid.electron_temp_k,
            background_grid.ion_temp_k,
            perturbed_density_m3,
            background_grid.neutral_species_cm3,
        )

    return IonosphereGrid(
        latitudes_deg=background_grid.latitudes_deg,
        longitudes_deg=background_grid.longitudes_deg,
        altitudes_km=background_grid.altitudes_km,
        iono_en_grid=perturbed_density_m3 / 1e6,
        iono_en_grid_5=perturbed_density_m3 / 1e6,
        collision_freq=collision_freq,
        iono_grid_parms=background_grid.iono_grid_parms,
        Bx=background_grid.Bx,
        By=background_grid.By,
        Bz=background_grid.Bz,
        geomag_grid_parms=background_grid.geomag_grid_parms,
        electron_temp_k=background_grid.electron_temp_k,
        ion_temp_k=background_grid.ion_temp_k,
        neutral_temp_k=background_grid.neutral_temp_k,
        neutral_species_cm3=background_grid.neutral_species_cm3,
        metadata=dict(background_grid.metadata),
    )


def _ray_path_ecef_m(ray: RayTrace) -> np.ndarray:
    heights = np.asarray(ray.path.get("height", []), dtype=np.float64)
    lats = np.asarray(ray.path.get("lat", []), dtype=np.float64)
    lons = np.asarray(ray.path.get("lon", []), dtype=np.float64)
    valid = np.isfinite(heights) & np.isfinite(lats) & np.isfinite(lons) & (heights < 1e40)
    if np.count_nonzero(valid) < 2:
        return np.zeros((0, 3), dtype=np.float64)
    return np.asarray(llh_to_ecef(lats[valid], lons[valid], heights[valid] * 1000.0), dtype=np.float64)


def _unit_vector(vector: np.ndarray) -> np.ndarray | None:
    norm = float(np.linalg.norm(vector))
    if not math.isfinite(norm) or norm <= 0.0:
        return None
    return np.asarray(vector, dtype=np.float64) / norm


def _plane_basis(normal: np.ndarray, rx_ecef_m: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    normal_hat = _unit_vector(normal)
    if normal_hat is None:
        return None
    up_guess = _unit_vector(rx_ecef_m)
    if up_guess is None:
        up_guess = np.array([0.0, 0.0, 1.0], dtype=np.float64)
    basis_x = np.cross(normal_hat, up_guess)
    if float(np.linalg.norm(basis_x)) < 1e-8:
        basis_x = np.cross(normal_hat, np.array([0.0, 0.0, 1.0], dtype=np.float64))
    basis_x = _unit_vector(basis_x)
    if basis_x is None:
        return None
    basis_y = _unit_vector(np.cross(normal_hat, basis_x))
    if basis_y is None:
        return None
    return basis_x, basis_y


def _receive_plane_coordinates(ray: RayTrace, rx: GeoPoint, basis_x: np.ndarray, basis_y: np.ndarray, normal_hat: np.ndarray) -> np.ndarray | None:
    distance = ray_point_distance(ray.path, rx)
    if distance.closest_point_ecef_m is None or not math.isfinite(distance.distance_m):
        return None
    projected = np.asarray(distance.closest_point_ecef_m, dtype=np.float64)
    rx_ecef_m = np.asarray(llh_to_ecef(rx.lat_deg, rx.lon_deg, rx.alt_km * 1000.0), dtype=np.float64)
    if rx_ecef_m.ndim > 1:
        rx_ecef_m = rx_ecef_m[0]
    projected = projected - np.dot(projected - rx_ecef_m, normal_hat) * normal_hat
    delta = projected - rx_ecef_m
    return np.array([float(np.dot(delta, basis_x)), float(np.dot(delta, basis_y))], dtype=np.float64)


def _estimate_ray_tube_metrics(
    tracer: PointToPointRayTracer,
    grid: IonosphereGrid,
    radar: CollmSkiymetRadar,
    target: GeoPoint,
    central_ray: RayTrace,
    ray_tube_model: RayTubeModel,
) -> RayTubeMetrics | None:
    if not central_ray.home or central_ray.launch_elevation_deg is None or central_ray.launch_bearing_deg is None:
        return None
    path_xyz = _ray_path_ecef_m(central_ray)
    distance = ray_point_distance(central_ray.path, target)
    if path_xyz.shape[0] < 2 or distance.segment_index is None:
        return None

    seg_idx = min(distance.segment_index, path_xyz.shape[0] - 2)
    tangent = path_xyz[seg_idx + 1] - path_xyz[seg_idx]
    normal_hat = _unit_vector(tangent)
    if normal_hat is None:
        return None
    rx_ecef_m = np.asarray(llh_to_ecef(target.lat_deg, target.lon_deg, target.alt_km * 1000.0), dtype=np.float64)
    if rx_ecef_m.ndim > 1:
        rx_ecef_m = rx_ecef_m[0]
    basis = _plane_basis(normal_hat, rx_ecef_m)
    if basis is None:
        return None
    basis_x, basis_y = basis

    delta_el = max(ray_tube_model.delta_elevation_deg, 1e-4)
    delta_az = max(ray_tube_model.delta_bearing_deg, 1e-4)
    offsets = {
        "ep": (central_ray.launch_elevation_deg + delta_el, central_ray.launch_bearing_deg),
        "em": (central_ray.launch_elevation_deg - delta_el, central_ray.launch_bearing_deg),
        "ap": (central_ray.launch_elevation_deg, central_ray.launch_bearing_deg + delta_az),
        "am": (central_ray.launch_elevation_deg, central_ray.launch_bearing_deg - delta_az),
    }
    coords: dict[str, np.ndarray] = {}
    for key, (elev_deg, bear_deg) in offsets.items():
        neighbor = tracer.trace_launch(
            tx=radar.tx,
            rx=target,
            frequency_mhz=radar.frequency_mhz,
            grid=grid,
            elevation_deg=elev_deg,
            bearing_deg=bear_deg,
            nhops=1,
        )
        coord = _receive_plane_coordinates(neighbor, target, basis_x, basis_y, normal_hat)
        if coord is None:
            return None
        coords[key] = coord

    dcoord_de = (coords["ep"] - coords["em"]) / (2.0 * math.radians(delta_el))
    dcoord_da = (coords["ap"] - coords["am"]) / (2.0 * math.radians(delta_az))
    area_m2_per_rad2 = abs(float(dcoord_de[0] * dcoord_da[1] - dcoord_de[1] * dcoord_da[0]))
    if not math.isfinite(area_m2_per_rad2) or area_m2_per_rad2 <= 0.0:
        return None

    reference_range_km = central_ray.group_range_to_rx_km
    if reference_range_km is None or reference_range_km <= 0.0:
        reference_range_km = central_ray.geometric_dist_to_rx_km
    reference_area_m2 = max(((reference_range_km or 1.0) * 1000.0) ** 2, 1.0)
    gain_db = -10.0 * math.log10(area_m2_per_rad2 / reference_area_m2)
    return RayTubeMetrics(area_m2_per_rad2=area_m2_per_rad2, gain_db=gain_db)


def _predict_physical_snr_db(
    model: SNRModel,
    slant_range_km: float,
    elevation_deg: float,
    absorption_db: float | None,
    ray_tube_metrics: RayTubeMetrics | None,
) -> float:
    if not math.isfinite(slant_range_km) or slant_range_km <= 0.0 or elevation_deg <= 0.0:
        return float("nan")
    reference_area_m2 = max((model.reference_range_km * 1000.0) ** 2, 1.0)
    area_m2 = reference_area_m2
    if ray_tube_metrics is not None:
        area_m2 = max(ray_tube_metrics.area_m2_per_rad2, 1.0)
    else:
        area_m2 = max((slant_range_km * 1000.0) ** 2, 1.0)
    path_term_db = -10.0 * math.log10(area_m2 / reference_area_m2)
    beam_ratio = max(
        math.sin(math.radians(elevation_deg)) / max(math.sin(math.radians(model.reference_elevation_deg)), 1e-6),
        1e-6,
    )
    beam_term_db = 10.0 * model.beam_power_exponent * math.log10(beam_ratio)
    absorption_term_db = -max(float(absorption_db or 0.0), 0.0)
    return max(model.floor_snr_db, model.reference_snr_db + path_term_db + beam_term_db + absorption_term_db)


def _solve_target_ray(
    radar: CollmSkiymetRadar,
    target: GeoPoint,
    *,
    homing_tolerance_m: float,
    tracer: PointToPointRayTracer | None = None,
    grid: IonosphereGrid | None = None,
    when_utc: dt.datetime | None = None,
    f107: float | None = None,
    d_region_model: str = "fpt2018",
    start_angles_deg: tuple[float, float] | None = None,
) -> RayTrace:
    if grid is None:
        if when_utc is None:
            raise ValueError("when_utc is required when solving a ray without a prebuilt grid")
        grid = build_pyiri_grid(
            when_utc.astimezone(dt.timezone.utc).replace(tzinfo=None),
            radar.tx,
            target,
            f107=f107,
            alt_min_km=60.0,
            alt_max_km=max(700.0, target.alt_km + 50.0),
            alt_step_km=5.0,
            lat_step_deg=0.5,
            lon_step_deg=0.5,
            lat_margin_deg=3.0,
            lon_margin_deg=3.0,
            d_region_model=d_region_model,
        )
    tracer = tracer or PointToPointRayTracer()
    if tracer is None:
        raise ValueError("tracer backend is unavailable")
    return tracer.trace_link(
        tx=radar.tx,
        rx=target,
        frequency_mhz=radar.frequency_mhz,
        grid=grid,
        nhops=1,
        homing_tolerance_m=homing_tolerance_m,
        start_angles_deg=start_angles_deg,
    )


def build_fof2_map_snapshot(
    series: PredictionSeries,
    *,
    lat_step_deg: float = 0.25,
    lon_step_deg: float = 0.25,
    margin_deg: float = 6.0,
    alt_min_km: float = 150.0,
    alt_max_km: float = 500.0,
    alt_step_km: float = 10.0,
) -> FoF2MapSnapshot:
    import PyIRI
    from PyIRI import main_library

    if not series.points:
        raise ValueError("cannot build a foF2 map without any predicted points")

    cpa_index = _closest_approach_index(series.points)
    cpa_point = series.points[cpa_index]
    snapshot_time = cpa_point.time_utc.astimezone(dt.timezone.utc).replace(tzinfo=None)
    elapsed_seconds = (
        (cpa_point.time_utc - series.config.start_utc).total_seconds()
        + series.config.ephemeris_time_shift_seconds
    )
    indices = resolve_space_weather_indices(
        snapshot_time,
        f107=series.config.f107,
        f107a=series.config.f107a,
    )

    track_lats = np.asarray([point.satellite.lat_deg for point in series.points], dtype=np.float64)
    track_lons_wrapped = np.asarray([wrap_longitude(point.satellite.lon_deg) for point in series.points], dtype=np.float64)
    reference_lon = float(np.mean(track_lons_wrapped)) if track_lons_wrapped.size else float(wrap_longitude(series.radar.tx.lon_deg))
    track_lons = _coerce_longitudes_near_reference(track_lons_wrapped, reference_lon)
    tx_lon = float(_coerce_longitudes_near_reference(np.array([wrap_longitude(series.radar.tx.lon_deg)], dtype=np.float64), reference_lon)[0])
    tx_ground = GeoPoint(series.radar.tx.lat_deg, series.radar.tx.lon_deg, 0.0)
    cpa_ground = GeoPoint(cpa_point.satellite.lat_deg, cpa_point.satellite.lon_deg, 0.0)
    track_bearing_deg = initial_bearing_deg(tx_ground, cpa_ground)
    wave_bearing_deg = _resolve_wave_bearing_deg(series.config.hooke_wave, track_bearing_deg)

    display_lat_min_deg = max(-90.0, float(np.min(np.append(track_lats, series.radar.tx.lat_deg))) - 1.0)
    display_lat_max_deg = min(90.0, float(np.max(np.append(track_lats, series.radar.tx.lat_deg))) + 1.0)
    display_lon_min_deg = float(np.min(track_lons))
    display_lon_max_deg = float(np.max(track_lons))
    if math.isclose(display_lon_min_deg, display_lon_max_deg, abs_tol=1e-9):
        display_lon_min_deg -= 0.5
        display_lon_max_deg += 0.5

    lat_min = max(-90.0, display_lat_min_deg - max(0.0, margin_deg - 1.0))
    lat_max = min(90.0, display_lat_max_deg + max(0.0, margin_deg - 1.0))
    lon_min = display_lon_min_deg
    lon_max = display_lon_max_deg

    lat_start = math.floor(lat_min / lat_step_deg) * lat_step_deg
    lat_stop = math.ceil(lat_max / lat_step_deg) * lat_step_deg
    lon_start = math.floor(lon_min / lon_step_deg) * lon_step_deg
    lon_stop = math.ceil(lon_max / lon_step_deg) * lon_step_deg

    latitudes_deg = np.arange(lat_start, lat_stop + 0.5 * lat_step_deg, lat_step_deg, dtype=np.float64)
    longitudes_deg = np.arange(lon_start, lon_stop + 0.5 * lon_step_deg, lon_step_deg, dtype=np.float64)
    altitudes_km = np.arange(alt_min_km, alt_max_km + 0.5 * alt_step_km, alt_step_km, dtype=np.float64)
    if latitudes_deg.size < 2 or longitudes_deg.size < 2 or altitudes_km.size < 2:
        raise ValueError("foF2 map axes must each contain at least two samples")

    lon_mesh, lat_mesh = np.meshgrid(longitudes_deg, latitudes_deg, indexing="xy")
    sample_lons = np.asarray([wrap_longitude(lon) for lon in lon_mesh.ravel()], dtype=np.float64)
    sample_lats = lat_mesh.ravel().astype(np.float64, copy=False)
    ut_hours = (
        snapshot_time.hour
        + snapshot_time.minute / 60.0
        + snapshot_time.second / 3600.0
        + snapshot_time.microsecond / 3.6e9
    )
    _, _, _, _, _, _, edens_m3 = main_library.IRI_density_1day(
        snapshot_time.year,
        snapshot_time.month,
        snapshot_time.day,
        np.array([ut_hours], dtype=np.float64),
        sample_lons,
        sample_lats,
        altitudes_km,
        indices.f107,
        PyIRI.coeff_dir,
    )
    density_alt_lat_lon_m3 = edens_m3[0].reshape(len(altitudes_km), len(latitudes_deg), len(longitudes_deg))
    peak_density_m3 = np.max(density_alt_lat_lon_m3, axis=0)
    background_fof2_mhz = 8.98e-6 * np.sqrt(np.maximum(peak_density_m3, 0.0))
    phase_coordinate_grid_km = _along_track_coordinate_km(
        tx_ground,
        wave_bearing_deg,
        lat_mesh,
        lon_mesh,
    )
    peak_density_m3 = _apply_hooke_wave_to_density(
        peak_density_m3[None, :, :],
        phase_coordinate_grid_km,
        np.array([series.config.hooke_wave.vertical_center_km], dtype=np.float64),
        elapsed_seconds,
        series.config.hooke_wave,
    )[0]
    fof2_mhz = 8.98e-6 * np.sqrt(np.maximum(peak_density_m3, 0.0))

    return FoF2MapSnapshot(
        time_utc=cpa_point.time_utc,
        latitudes_deg=latitudes_deg,
        longitudes_deg=longitudes_deg,
        fof2_mhz=fof2_mhz.astype(np.float32, copy=False),
        background_fof2_mhz=background_fof2_mhz.astype(np.float32, copy=False),
        wave_origin_lat_deg=tx_ground.lat_deg,
        wave_origin_lon_deg=tx_ground.lon_deg,
        wave_bearing_deg=wave_bearing_deg,
        wave_elapsed_seconds=elapsed_seconds,
        hooke_wave=series.config.hooke_wave,
        track_latitudes_deg=track_lats,
        track_longitudes_deg=track_lons,
        transmitter_lat_deg=series.radar.tx.lat_deg,
        transmitter_lon_deg=tx_lon,
        cpa_lat_deg=cpa_point.satellite.lat_deg,
        cpa_lon_deg=float(track_lons[cpa_index]),
        cpa_index=cpa_index,
        f107=indices.f107,
        f107a=indices.f107a,
        display_lon_min_deg=display_lon_min_deg,
        display_lon_max_deg=display_lon_max_deg,
        display_lat_min_deg=display_lat_min_deg,
        display_lat_max_deg=display_lat_max_deg,
    )


def _ray_subpoint_range_from_tx_km(
    tx_ground: GeoPoint,
    path_lats_deg: np.ndarray,
    path_lons_deg: np.ndarray,
) -> np.ndarray:
    ray_lats = np.asarray(path_lats_deg, dtype=np.float64)
    ray_lons = np.asarray(path_lons_deg, dtype=np.float64)
    out = np.empty(ray_lats.shape, dtype=np.float64)
    for idx, (lat_deg, lon_deg) in enumerate(zip(ray_lats, ray_lons)):
        out[idx] = _surface_distance_km(tx_ground, GeoPoint(float(lat_deg), float(lon_deg), 0.0))
    return np.maximum.accumulate(out)


def build_cpa_swath_snapshot(
    series: PredictionSeries,
    *,
    waypoint_count: int = 256,
    alt_min_km: float = 0.0,
    alt_max_km: float = 500.0,
    alt_step_km: float = 5.0,
) -> SwathSnapshot:
    import PyIRI
    from PyIRI import main_library

    if not series.points:
        raise ValueError("cannot build a swath without any predicted points")

    cpa_index = _closest_approach_index(series.points)
    cpa_point = series.points[cpa_index]
    snapshot_time = cpa_point.time_utc.astimezone(dt.timezone.utc).replace(tzinfo=None)
    elapsed_seconds = (
        (cpa_point.time_utc - series.config.start_utc).total_seconds()
        + series.config.ephemeris_time_shift_seconds
    )
    indices = resolve_space_weather_indices(
        snapshot_time,
        f107=series.config.f107,
        f107a=series.config.f107a,
    )

    tx_ground = GeoPoint(series.radar.tx.lat_deg, series.radar.tx.lon_deg, 0.0)
    rx_ground = GeoPoint(cpa_point.satellite.lat_deg, cpa_point.satellite.lon_deg, 0.0)
    track_lats_deg, track_lons_deg = great_circle_waypoints(tx_ground, rx_ground, count=waypoint_count)
    surface_range_km = _surface_distance_km(tx_ground, rx_ground)
    along_track_km = np.linspace(0.0, surface_range_km, waypoint_count, dtype=np.float64)
    wave_bearing_deg = _resolve_wave_bearing_deg(series.config.hooke_wave, initial_bearing_deg(tx_ground, rx_ground))
    phase_coordinate_km = _along_track_coordinate_km(
        tx_ground,
        wave_bearing_deg,
        np.asarray(track_lats_deg, dtype=np.float64),
        np.asarray(track_lons_deg, dtype=np.float64),
    )

    altitudes_km = np.arange(alt_min_km, alt_max_km + 0.5 * alt_step_km, alt_step_km, dtype=np.float64)
    if altitudes_km.size < 2:
        raise ValueError("swath altitude axis must contain at least two samples")

    wrapped_track_lons_deg = np.asarray([wrap_longitude(lon) for lon in track_lons_deg], dtype=np.float64)
    ut_hours = (
        snapshot_time.hour
        + snapshot_time.minute / 60.0
        + snapshot_time.second / 3600.0
        + snapshot_time.microsecond / 3.6e9
    )
    _, _, _, _, _, _, edens_m3 = main_library.IRI_density_1day(
        snapshot_time.year,
        snapshot_time.month,
        snapshot_time.day,
        np.array([ut_hours], dtype=np.float64),
        wrapped_track_lons_deg,
        np.asarray(track_lats_deg, dtype=np.float64),
        altitudes_km,
        indices.f107,
        PyIRI.coeff_dir,
    )
    density_alt_track_m3 = edens_m3[0].reshape(len(altitudes_km), waypoint_count)
    density_alt_track_m3 = _apply_hooke_wave_to_density(
        density_alt_track_m3,
        phase_coordinate_km,
        altitudes_km,
        elapsed_seconds,
        series.config.hooke_wave,
    )
    plasma_freq_mhz = 8.98e-6 * np.sqrt(np.maximum(density_alt_track_m3, 0.0))

    cpa_context = _build_real_pass_context(
        series.radar,
        [cpa_point.satellite],
        cpa_point.time_utc,
        series.config,
        wave_bearing_deg,
    )
    cpa_grid = _wave_perturbed_grid(
        cpa_context.background_grid,
        cpa_context.phase_coordinate_km,
        elapsed_seconds,
        series.config.hooke_wave,
    )
    ray = _solve_target_ray(
        series.radar,
        cpa_point.satellite,
        homing_tolerance_m=series.config.homing_tolerance_m,
        tracer=cpa_context.tracer,
        grid=cpa_grid,
    )
    ray_lats_deg = np.asarray(ray.path.get("lat", []), dtype=np.float64)
    ray_lons_deg = np.asarray(ray.path.get("lon", []), dtype=np.float64)
    ray_altitudes_km = np.asarray(ray.path.get("height", []), dtype=np.float64)
    valid = (
        np.isfinite(ray_lats_deg)
        & np.isfinite(ray_lons_deg)
        & np.isfinite(ray_altitudes_km)
        & (ray_altitudes_km < 1e40)
    )
    ray_lats_deg = ray_lats_deg[valid]
    ray_lons_deg = ray_lons_deg[valid]
    ray_altitudes_km = ray_altitudes_km[valid]
    ray_along_track_km = _ray_subpoint_range_from_tx_km(
        tx_ground,
        ray_lats_deg,
        ray_lons_deg,
    )

    return SwathSnapshot(
        time_utc=cpa_point.time_utc,
        along_track_km=along_track_km,
        altitudes_km=altitudes_km,
        plasma_freq_mhz=plasma_freq_mhz.astype(np.float32, copy=False),
        ray_along_track_km=ray_along_track_km.astype(np.float32, copy=False),
        ray_altitudes_km=ray_altitudes_km.astype(np.float32, copy=False),
        tx_ground_label=series.radar.site,
        rx_ground_label="CPA subpoint",
        surface_range_km=surface_range_km,
        target_alt_km=cpa_point.satellite.alt_km,
        f107=indices.f107,
        f107a=indices.f107a,
    )


def _build_point_prediction(
    radar: CollmSkiymetRadar,
    when: dt.datetime,
    start_utc: dt.datetime,
    target: GeoPoint,
    absolute_range_km: float,
    folded_range_km: float,
    absolute_doppler_hz: float,
    aliased_doppler_hz: float,
    homing_tolerance_m: float,
    snr_model: SNRModel,
    ephemeris_time_shift_seconds: float,
    hooke_wave: HookeWaveModel,
    wave_bearing_deg: float,
    ray_tube_model: RayTubeModel,
    real_context: RealPassContext | None = None,
    start_angles_deg: tuple[float, float] | None = None,
    f107: float | None = None,
    real_d_region_model: str = "fpt2018",
) -> PredictionPoint:
    los_bearing_deg, los_elevation_deg, slant_range_km = _line_of_sight_angles(radar.tx, target)
    tx_ground = GeoPoint(radar.tx.lat_deg, radar.tx.lon_deg, 0.0)
    phase_coordinate_km = float(_along_track_coordinate_km(tx_ground, wave_bearing_deg, target.lat_deg, target.lon_deg))
    elapsed_seconds = (when - start_utc).total_seconds() + ephemeris_time_shift_seconds

    launch_bearing_deg: float | None = None
    launch_elevation_deg: float | None = None
    group_range_km: float | None = None
    total_absorption_db: float | None = None
    ray_tube_area_m2: float | None = None
    ray_tube_gain_db: float | None = None
    max_path_alt_km: float | None = None
    home = False
    error_m: float | None = None

    if los_elevation_deg >= 0.0:
        if real_context is not None:
            sample_grid = _wave_perturbed_grid(
                real_context.background_grid,
                real_context.phase_coordinate_km,
                elapsed_seconds,
                hooke_wave,
            )
            ray = _solve_target_ray(
                radar,
                target,
                homing_tolerance_m=homing_tolerance_m,
                tracer=real_context.tracer,
                grid=sample_grid,
                start_angles_deg=start_angles_deg,
            )
            ray_tube_metrics = _estimate_ray_tube_metrics(
                real_context.tracer,
                sample_grid,
                radar,
                target,
                ray,
                ray_tube_model,
            )
        else:
            ray = _solve_target_ray(
                radar,
                target,
                homing_tolerance_m=homing_tolerance_m,
                when_utc=when,
                f107=f107,
                d_region_model=real_d_region_model,
                start_angles_deg=start_angles_deg,
            )
            ray_tube_metrics = None
        launch_bearing_deg = ray.launch_bearing_deg
        launch_elevation_deg = ray.launch_elevation_deg
        group_range_km = ray.group_range_to_rx_km
        total_absorption_db = ray.total_absorption_db
        error_m = ray.error_m if math.isfinite(ray.error_m) else None
        home = ray.home
        if ray_tube_metrics is not None:
            ray_tube_area_m2 = ray_tube_metrics.area_m2_per_rad2
            ray_tube_gain_db = ray_tube_metrics.gain_db
        heights = np.asarray(ray.path.get("height", []), dtype=float)
        heights = heights[np.isfinite(heights) & (heights < 1e40)]
        if heights.size:
            max_path_alt_km = float(np.max(heights))

    predicted_peak_snr_db = _predict_physical_snr_db(
        snr_model,
        slant_range_km=slant_range_km,
        elevation_deg=los_elevation_deg,
        absorption_db=total_absorption_db,
        ray_tube_metrics=None if ray_tube_area_m2 is None else RayTubeMetrics(ray_tube_area_m2, ray_tube_gain_db or 0.0),
    )
    return PredictionPoint(
        time_utc=when,
        satellite=target,
        los_bearing_deg=los_bearing_deg,
        los_elevation_deg=los_elevation_deg,
        absolute_range_km=absolute_range_km,
        folded_range_km=folded_range_km,
        absolute_doppler_hz=absolute_doppler_hz,
        aliased_doppler_hz=aliased_doppler_hz,
        predicted_peak_snr_db=predicted_peak_snr_db,
        launch_bearing_deg=launch_bearing_deg,
        launch_elevation_deg=launch_elevation_deg,
        group_range_km=group_range_km,
        total_absorption_db=total_absorption_db,
        ray_tube_area_m2=ray_tube_area_m2,
        ray_tube_gain_db=ray_tube_gain_db,
        max_path_alt_km=max_path_alt_km,
        home=home,
        error_m=error_m,
    )


def _build_synthetic_delay_image(
    radar: CollmSkiymetRadar,
    points: Sequence[PredictionPoint],
    noise_seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    bin_size_km = C_KM_PER_S / radar.channel_rate_hz
    range_axis_km = np.arange(0.0, radar.unambiguous_range_km + 0.5 * bin_size_km, bin_size_km, dtype=np.float64)
    if range_axis_km.size < 2:
        range_axis_km = np.array([0.0, radar.unambiguous_range_km], dtype=np.float64)

    rng = np.random.default_rng(noise_seed)
    image = rng.normal(loc=-0.5, scale=0.12, size=(range_axis_km.size, len(points))).astype(np.float32)
    sigma_km = 1.4 * bin_size_km

    for col, point in enumerate(points):
        if not math.isfinite(point.folded_range_km) or not math.isfinite(point.predicted_peak_snr_db):
            continue
        center = point.folded_range_km
        ridge = point.predicted_peak_snr_db - 2.8 * ((range_axis_km - center) / max(sigma_km, 1e-6)) ** 2
        image[:, col] = np.maximum(image[:, col], ridge.astype(np.float32))
    return range_axis_km, image


def predict_collm_skiymet_time_series(
    config: TimeSeriesConfig,
    radar: CollmSkiymetRadar | None = None,
) -> PredictionSeries:
    radar = radar or CollmSkiymetRadar()
    sat = load_tle(config.tle_file.expanduser())
    row_times = _time_grid(config.start_utc.astimezone(dt.timezone.utc), config.seconds, config.step_seconds)
    unix_times = np.asarray([time.timestamp() for time in row_times], dtype=np.float64)
    shifted_unix_times = unix_times + config.ephemeris_time_shift_seconds

    radar_ecef_m = np.asarray(
        llh_to_ecef(radar.tx.lat_deg, radar.tx.lon_deg, radar.tx.alt_km * 1000.0),
        dtype=np.float64,
    )
    if radar_ecef_m.ndim > 1:
        radar_ecef_m = radar_ecef_m[0]

    delays_s, absolute_dopplers_hz = predict_delay_doppler(sat, shifted_unix_times, radar_ecef_m, radar.freq_hz)
    absolute_ranges_km = delays_s * C_KM_PER_S
    folded_ranges_km = np.mod(absolute_ranges_km + config.range_offset_km, radar.unambiguous_range_km)
    aliased_dopplers_hz = _alias_to_prf_band(absolute_dopplers_hz + config.doppler_offset_hz, radar.prf_hz)
    targets = [_ecef_to_geo(_sat_ecef_at_unix(sat, float(unix_time_s))) for unix_time_s in shifted_unix_times]
    finite_range_indices = np.flatnonzero(np.isfinite(absolute_ranges_km))
    cpa_index = int(finite_range_indices[np.argmin(absolute_ranges_km[finite_range_indices])]) if finite_range_indices.size else 0
    tx_ground = GeoPoint(radar.tx.lat_deg, radar.tx.lon_deg, 0.0)
    cpa_ground = GeoPoint(targets[cpa_index].lat_deg, targets[cpa_index].lon_deg, 0.0)
    wave_bearing_deg = _resolve_wave_bearing_deg(config.hooke_wave, initial_bearing_deg(tx_ground, cpa_ground))
    real_context = _build_real_pass_context(
        radar,
        targets,
        row_times[cpa_index],
        config,
        wave_bearing_deg,
    )

    points: list[PredictionPoint] = []
    start_angles_deg: tuple[float, float] | None = None
    for when, abs_range_km, folded_range_km, abs_doppler_hz, aliased_doppler_hz, target in zip(
        row_times,
        absolute_ranges_km,
        folded_ranges_km,
        absolute_dopplers_hz,
        aliased_dopplers_hz,
        targets,
    ):
        point = _build_point_prediction(
            radar=radar,
            when=when,
            start_utc=config.start_utc,
            target=target,
            absolute_range_km=float(abs_range_km),
            folded_range_km=float(folded_range_km),
            absolute_doppler_hz=float(abs_doppler_hz),
            aliased_doppler_hz=float(aliased_doppler_hz),
            homing_tolerance_m=config.homing_tolerance_m,
            snr_model=config.snr_model,
            ephemeris_time_shift_seconds=config.ephemeris_time_shift_seconds,
            hooke_wave=config.hooke_wave,
            wave_bearing_deg=wave_bearing_deg,
            ray_tube_model=config.ray_tube_model,
            real_context=real_context,
            start_angles_deg=start_angles_deg,
            f107=config.f107,
            real_d_region_model=config.real_d_region_model,
        )
        if point.launch_elevation_deg is not None and point.launch_bearing_deg is not None:
            start_angles_deg = (point.launch_elevation_deg, point.launch_bearing_deg)
        if point.los_elevation_deg < config.min_elevation_deg:
            point = PredictionPoint(
                time_utc=point.time_utc,
                satellite=point.satellite,
                los_bearing_deg=point.los_bearing_deg,
                los_elevation_deg=point.los_elevation_deg,
                absolute_range_km=point.absolute_range_km,
                folded_range_km=point.folded_range_km,
                absolute_doppler_hz=point.absolute_doppler_hz,
                aliased_doppler_hz=point.aliased_doppler_hz,
                predicted_peak_snr_db=float("nan"),
                launch_bearing_deg=None,
                launch_elevation_deg=None,
                group_range_km=None,
                total_absorption_db=None,
                ray_tube_area_m2=None,
                ray_tube_gain_db=None,
                max_path_alt_km=None,
                home=False,
                error_m=None,
            )
        points.append(point)

    range_axis_km, image_snr_db = _build_synthetic_delay_image(radar, points, config.noise_seed)
    return PredictionSeries(
        radar=radar,
        config=config,
        points=tuple(points),
        range_axis_km=range_axis_km,
        image_snr_db=image_snr_db,
    )


def write_prediction_csv(path: Path, series: PredictionSeries) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "utc_time",
                "sat_lat_deg",
                "sat_lon_deg",
                "sat_alt_km",
                "los_bearing_deg",
                "los_elevation_deg",
                "absolute_range_km",
                "folded_range_km",
                "absolute_doppler_hz",
                "aliased_doppler_hz",
                "predicted_peak_snr_db",
                "launch_bearing_deg",
                "launch_elevation_deg",
                "group_range_km",
                "total_absorption_db",
                "ray_tube_area_m2",
                "ray_tube_gain_db",
                "max_path_alt_km",
                "home",
                "error_m",
            ]
        )
        for point in series.points:
            writer.writerow(
                [
                    point.time_utc.isoformat(),
                    f"{point.satellite.lat_deg:.6f}",
                    f"{wrap_longitude(point.satellite.lon_deg):.6f}",
                    f"{point.satellite.alt_km:.3f}",
                    f"{point.los_bearing_deg:.3f}",
                    f"{point.los_elevation_deg:.3f}",
                    f"{point.absolute_range_km:.3f}",
                    f"{point.folded_range_km:.3f}",
                    f"{point.absolute_doppler_hz:.3f}",
                    f"{point.aliased_doppler_hz:.3f}",
                    "" if not math.isfinite(point.predicted_peak_snr_db) else f"{point.predicted_peak_snr_db:.3f}",
                    "" if point.launch_bearing_deg is None else f"{point.launch_bearing_deg:.3f}",
                    "" if point.launch_elevation_deg is None else f"{point.launch_elevation_deg:.3f}",
                    "" if point.group_range_km is None else f"{point.group_range_km:.3f}",
                    "" if point.total_absorption_db is None else f"{point.total_absorption_db:.3f}",
                    "" if point.ray_tube_area_m2 is None else f"{point.ray_tube_area_m2:.6e}",
                    "" if point.ray_tube_gain_db is None else f"{point.ray_tube_gain_db:.3f}",
                    "" if point.max_path_alt_km is None else f"{point.max_path_alt_km:.3f}",
                    int(point.home),
                    "" if point.error_m is None else f"{point.error_m:.3f}",
                ]
            )


def plot_prediction(
    path: Path,
    series: PredictionSeries,
    *,
    measurements: Sequence[MeasurementPoint] = (),
) -> Path:
    import matplotlib.dates as mdates  # type: ignore
    import matplotlib.pyplot as plt  # type: ignore
    from matplotlib.gridspec import GridSpec  # type: ignore

    radar = series.radar
    times = [point.time_utc for point in series.points]
    x_nums = mdates.date2num(times)
    if len(x_nums) > 1:
        step = float(np.median(np.diff(x_nums)))
    else:
        step = 1.0 / 86400.0
    x_edges = np.empty(len(x_nums) + 1, dtype=np.float64)
    if len(x_nums) == 1:
        x_edges[0] = float(x_nums[0] - 0.5 * step)
        x_edges[1] = float(x_nums[0] + 0.5 * step)
    else:
        x_edges[1:-1] = 0.5 * (x_nums[:-1] + x_nums[1:])
        x_edges[0] = float(x_nums[0] - 0.5 * (x_nums[1] - x_nums[0]))
        x_edges[-1] = float(x_nums[-1] + 0.5 * (x_nums[-1] - x_nums[-2]))
    range_edges = np.empty(series.range_axis_km.size + 1, dtype=np.float64)
    if series.range_axis_km.size == 1:
        bin_width = C_KM_PER_S / radar.channel_rate_hz
        range_edges[0] = series.range_axis_km[0] - 0.5 * bin_width
        range_edges[1] = series.range_axis_km[0] + 0.5 * bin_width
    else:
        range_edges[1:-1] = 0.5 * (series.range_axis_km[:-1] + series.range_axis_km[1:])
        range_edges[0] = series.range_axis_km[0] - 0.5 * (series.range_axis_km[1] - series.range_axis_km[0])
        range_edges[-1] = series.range_axis_km[-1] + 0.5 * (series.range_axis_km[-1] - series.range_axis_km[-2])

    aliased_dopplers = np.asarray([point.aliased_doppler_hz for point in series.points], dtype=np.float64)
    peak_snr_db = np.asarray([point.predicted_peak_snr_db for point in series.points], dtype=np.float64)
    folded_ranges = np.asarray([point.folded_range_km for point in series.points], dtype=np.float64)
    measurement_times = [point.time_utc for point in measurements]
    measured_ranges = np.asarray([point.peak_range_km for point in measurements], dtype=np.float64)
    measured_dopplers = np.asarray([point.doppler_hz for point in measurements], dtype=np.float64)
    measured_snrs = np.asarray([point.peak_snr_db for point in measurements], dtype=np.float64)

    fig = plt.figure(figsize=(13.5, 11.8))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[4.0, 1.1, 1.0], width_ratios=[40.0, 1.6], hspace=0.18, wspace=0.08)
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[1, 0], sharex=ax0)
    ax2 = fig.add_subplot(gs[2, 0], sharex=ax0)
    cax = fig.add_subplot(gs[:, 1])

    def _apply_grid(ax) -> None:
        ax.set_axisbelow(True)
        ax.minorticks_on()
        ax.grid(True, which="major", alpha=0.28, linewidth=0.7)
        ax.grid(True, which="minor", alpha=0.12, linewidth=0.5, linestyle=":")

    cf = ax0.pcolormesh(
        x_edges,
        range_edges,
        np.ma.masked_invalid(series.image_snr_db),
        cmap="viridis",
        shading="auto",
        vmin=-2.5,
        vmax=max(9.0, float(np.nanmax(series.image_snr_db)) if series.image_snr_db.size else 9.0),
    )
    ax0.plot(times, folded_ranges, color="white", linewidth=2.0, alpha=0.9)
    if measurement_times:
        ax0.scatter(measurement_times, measured_ranges, s=22, color="tab:red", edgecolors="white", linewidths=0.4, zorder=5)
    ax0.set_ylabel("One-way delay-equivalent range (km)", fontsize=15)
    ax0.set_ylim(0.0, radar.unambiguous_range_km)
    ax0.set_title(radar.title, fontsize=21, fontweight="bold")
    ax0.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    ax0.tick_params(axis="x", labelbottom=False)
    _apply_grid(ax0)

    finite_doppler = np.isfinite(aliased_dopplers)
    if np.any(finite_doppler):
        ax1.plot(np.asarray(times, dtype=object)[finite_doppler], aliased_dopplers[finite_doppler], color="tab:orange", linewidth=1.6)
        ax1.scatter(np.asarray(times, dtype=object)[finite_doppler], aliased_dopplers[finite_doppler], s=14, color="tab:orange", linewidths=0.0)
    if measurement_times:
        ax1.scatter(measurement_times, measured_dopplers, s=20, color="tab:red", edgecolors="white", linewidths=0.4, zorder=5)
    ax1.axhline(0.0, color="white", linewidth=0.8, alpha=0.5)
    ax1.set_ylabel("Doppler (Hz)", fontsize=13)
    ax1.set_ylim(-radar.doppler_nyquist_hz, radar.doppler_nyquist_hz)
    ax1.tick_params(axis="x", labelbottom=False)
    _apply_grid(ax1)

    finite_peak = np.isfinite(peak_snr_db)
    if np.any(finite_peak):
        ax2.plot(np.asarray(times, dtype=object)[finite_peak], peak_snr_db[finite_peak], color="tab:blue", linewidth=1.4, alpha=0.95)
        ax2.scatter(np.asarray(times, dtype=object)[finite_peak], peak_snr_db[finite_peak], s=10, color="tab:blue", alpha=0.8, linewidths=0.0)
        if measurement_times:
            ax2.scatter(measurement_times, measured_snrs, s=20, color="tab:red", edgecolors="white", linewidths=0.4, zorder=5)
        lo, hi = np.percentile(peak_snr_db[finite_peak], [2.0, 98.0])
        if not math.isfinite(lo) or not math.isfinite(hi) or lo == hi:
            lo = float(np.nanmin(peak_snr_db[finite_peak]))
            hi = float(np.nanmax(peak_snr_db[finite_peak]))
        if lo == hi:
            lo -= 1.0
            hi += 1.0
        pad = max(0.1 * (hi - lo), 1.0)
        ax2.set_ylim(lo - pad, hi + pad)
    ax2.set_ylabel("Peak SNR (dB)", fontsize=13)
    ax2.set_xlabel("UTC time", fontsize=15)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
    ax2.tick_params(axis="x", labelrotation=20)
    _apply_grid(ax2)

    cbar = fig.colorbar(cf, cax=cax, label="Matched-filter power SNR (dB)")
    cbar.ax.tick_params(labelsize=12)

    fig.autofmt_xdate()
    fig.subplots_adjust(left=0.08, right=0.95, top=0.92, bottom=0.10)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)
    return path


def plot_fof2_map(
    path: Path,
    snapshot: FoF2MapSnapshot,
    radar: CollmSkiymetRadar,
    *,
    mode: str = "total",
) -> Path:
    from scipy.interpolate import RegularGridInterpolator
    import matplotlib.pyplot as plt  # type: ignore

    lat_step = float(np.median(np.diff(snapshot.latitudes_deg))) if snapshot.latitudes_deg.size > 1 else 0.25
    lon_step = float(np.median(np.diff(snapshot.longitudes_deg))) if snapshot.longitudes_deg.size > 1 else 0.25
    interp_factor = 8.0
    dense_lat_step = lat_step / interp_factor
    dense_lon_step = lon_step / interp_factor

    dense_latitudes_deg = np.arange(
        snapshot.display_lat_min_deg,
        snapshot.display_lat_max_deg + 0.5 * dense_lat_step,
        dense_lat_step,
        dtype=np.float64,
    )
    dense_longitudes_deg = np.arange(
        snapshot.display_lon_min_deg,
        snapshot.display_lon_max_deg + 0.5 * dense_lon_step,
        dense_lon_step,
        dtype=np.float64,
    )
    dense_latitudes_deg = np.clip(dense_latitudes_deg, snapshot.latitudes_deg[0], snapshot.latitudes_deg[-1])
    dense_longitudes_deg = np.clip(dense_longitudes_deg, snapshot.longitudes_deg[0], snapshot.longitudes_deg[-1])
    dense_latitudes_deg = np.unique(dense_latitudes_deg)
    dense_longitudes_deg = np.unique(dense_longitudes_deg)
    dense_lon_centers, dense_lat_centers = np.meshgrid(
        dense_longitudes_deg,
        dense_latitudes_deg,
        indexing="xy",
    )

    background_interpolator = RegularGridInterpolator(
        (snapshot.latitudes_deg, snapshot.longitudes_deg),
        snapshot.background_fof2_mhz,
        bounds_error=False,
        fill_value=np.nan,
    )
    query_points = np.column_stack((dense_lat_centers.ravel(), dense_lon_centers.ravel()))
    dense_background_fof2_mhz = background_interpolator(query_points).reshape(dense_lat_centers.shape)

    wave_origin = GeoPoint(snapshot.wave_origin_lat_deg, snapshot.wave_origin_lon_deg, 0.0)
    dense_phase_coordinate_km = _along_track_coordinate_km(
        wave_origin,
        snapshot.wave_bearing_deg,
        dense_lat_centers,
        dense_lon_centers,
    )
    dense_wave_fraction = _hooke_wave_fraction(
        snapshot.hooke_wave,
        dense_phase_coordinate_km,
        snapshot.hooke_wave.vertical_center_km,
        snapshot.wave_elapsed_seconds,
    )
    dense_fof2_ratio = np.sqrt(np.maximum(1.0 + dense_wave_fraction, 1e-6))

    if mode == "total":
        dense_fof2_mhz = dense_background_fof2_mhz * dense_fof2_ratio
        cmap = "viridis"
        colorbar_label = "foF2 (MHz)"
        title_suffix = "foF2"
    elif mode == "delta":
        dense_fof2_mhz = dense_background_fof2_mhz * (dense_fof2_ratio - 1.0)
        cmap = "RdBu_r"
        colorbar_label = "Delta foF2 (MHz)"
        title_suffix = "Delta foF2"
    elif mode == "fractional":
        dense_fof2_mhz = 100.0 * (dense_fof2_ratio - 1.0)
        cmap = "RdBu_r"
        colorbar_label = "Delta foF2 / background (%)"
        title_suffix = "Fractional foF2 perturbation"
    else:
        raise ValueError(f"unsupported foF2 map mode {mode!r}")

    lon_edges = _centers_to_edges(dense_longitudes_deg)
    lat_edges = _centers_to_edges(dense_latitudes_deg)
    lon_mesh, lat_mesh = np.meshgrid(lon_edges, lat_edges, indexing="xy")

    finite = dense_fof2_mhz[np.isfinite(dense_fof2_mhz)]
    if finite.size:
        vmin = float(np.min(finite))
        vmax = float(np.max(finite))
        if vmax <= vmin:
            pad = max(abs(vmin) * 1e-3, 1e-3)
            vmin -= pad
            vmax += pad
    else:
        vmin = 0.0
        vmax = 10.0

    fig, ax = plt.subplots(figsize=(11.5, 7.8))
    surface = ax.pcolormesh(
        lon_mesh,
        lat_mesh,
        dense_fof2_mhz,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="flat",
        antialiased=False,
        rasterized=True,
    )
    ax.plot(
        snapshot.track_longitudes_deg,
        snapshot.track_latitudes_deg,
        color="white",
        linewidth=2.6,
        alpha=0.95,
        label="ISS ground track",
    )
    ax.scatter(
        snapshot.track_longitudes_deg[0],
        snapshot.track_latitudes_deg[0],
        color="white",
        edgecolors="black",
        linewidths=0.6,
        s=42,
        zorder=4,
        label="Track start",
    )
    ax.scatter(
        snapshot.cpa_lon_deg,
        snapshot.cpa_lat_deg,
        color="tab:red",
        edgecolors="white",
        linewidths=0.8,
        s=68,
        zorder=5,
        label="CPA",
    )
    ax.scatter(
        snapshot.transmitter_lon_deg,
        snapshot.transmitter_lat_deg,
        marker="*",
        color="gold",
        edgecolors="black",
        linewidths=0.9,
        s=220,
        zorder=6,
        label=radar.site,
    )
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")
    ax.set_xlim(snapshot.display_lon_min_deg, snapshot.display_lon_max_deg)
    ax.set_ylim(snapshot.display_lat_min_deg, snapshot.display_lat_max_deg)
    ax.set_title(
        f"{radar.label} {title_suffix} at CPA ({snapshot.time_utc.isoformat()})",
        fontsize=16,
        fontweight="bold",
    )
    ax.grid(True, which="major", alpha=0.25, linewidth=0.7)
    ax.minorticks_on()
    ax.grid(True, which="minor", alpha=0.12, linewidth=0.5, linestyle=":")
    ax.legend(loc="best", framealpha=0.9)

    cbar = fig.colorbar(surface, ax=ax, pad=0.02, label=colorbar_label)
    cbar.ax.tick_params(labelsize=10)
    subtitle = (
        f"Track span: {snapshot.track_latitudes_deg[0]:.2f}° to {snapshot.track_latitudes_deg[-1]:.2f}° lat, "
        f"F10.7={snapshot.f107:.1f}, F10.7a={snapshot.f107a:.1f}"
    )
    ax.text(0.01, 0.01, subtitle, transform=ax.transAxes, ha="left", va="bottom", fontsize=9, color="white",
            bbox={"facecolor": "black", "alpha": 0.28, "pad": 4, "edgecolor": "none"})

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    return path


def plot_cpa_swath(path: Path, snapshot: SwathSnapshot, radar: CollmSkiymetRadar) -> Path:
    import matplotlib.pyplot as plt  # type: ignore

    along_edges = _centers_to_edges(snapshot.along_track_km)
    alt_edges = _centers_to_edges(snapshot.altitudes_km)
    along_mesh, alt_mesh = np.meshgrid(along_edges, alt_edges, indexing="xy")

    finite = snapshot.plasma_freq_mhz[np.isfinite(snapshot.plasma_freq_mhz)]
    if finite.size:
        vmin = float(np.min(finite))
        vmax = float(np.max(finite))
        if vmax <= vmin:
            pad = max(abs(vmin) * 1e-3, 1e-3)
            vmin -= pad
            vmax += pad
    else:
        vmin = 0.0
        vmax = 10.0

    fig, ax = plt.subplots(figsize=(12.0, 6.8))
    surface = ax.pcolormesh(
        along_mesh,
        alt_mesh,
        snapshot.plasma_freq_mhz,
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
        shading="flat",
        antialiased=False,
        rasterized=True,
    )
    ax.plot(
        snapshot.ray_along_track_km,
        snapshot.ray_altitudes_km,
        color="white",
        linewidth=2.8,
        alpha=0.95,
        label="Homed ray path",
    )
    ax.scatter(
        [0.0, snapshot.surface_range_km],
        [radar.tx.alt_km, snapshot.target_alt_km],
        s=56,
        c=["gold", "tab:red"],
        edgecolors="black",
        linewidths=0.7,
        zorder=5,
    )
    ax.set_xlim(0.0, snapshot.surface_range_km)
    ax.set_ylim(float(snapshot.altitudes_km[0]), float(snapshot.altitudes_km[-1]))
    ax.text(
        0.01 * snapshot.surface_range_km,
        float(snapshot.altitudes_km[0]) + 8.0,
        snapshot.tx_ground_label,
        ha="left",
        va="bottom",
        fontsize=10,
        color="white",
    )
    ax.text(
        0.985 * snapshot.surface_range_km,
        min(snapshot.target_alt_km + 8.0, float(snapshot.altitudes_km[-1]) - 10.0),
        "ISS at CPA",
        ha="right",
        va="bottom",
        fontsize=10,
        color="white",
    )
    ax.set_xlabel(f"Along-track distance from {snapshot.tx_ground_label} (km)")
    ax.set_ylabel("Altitude (km)")
    ax.set_title(
        f"{radar.label} CPA Swath and Ray Path ({snapshot.time_utc.isoformat()})",
        fontsize=16,
        fontweight="bold",
    )
    ax.set_axisbelow(True)
    ax.minorticks_on()
    ax.grid(True, which="major", alpha=0.28, linewidth=0.7)
    ax.grid(True, which="minor", alpha=0.12, linewidth=0.5, linestyle=":")
    ax.legend(loc="upper right", framealpha=0.9)

    cbar = fig.colorbar(surface, ax=ax, pad=0.02, label="Plasma frequency (MHz)")
    cbar.ax.tick_params(labelsize=10)
    ax.text(
        0.01,
        0.01,
        f"Ground track: {snapshot.surface_range_km:.1f} km, F10.7={snapshot.f107:.1f}, F10.7a={snapshot.f107a:.1f}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9,
        color="white",
        bbox={"facecolor": "black", "alpha": 0.28, "pad": 4, "edgecolor": "none"},
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)
    return path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Predict a Collm Skiymet-style range/Doppler/SNR time series from an ISS TLE.",
    )
    parser.add_argument("--tle-file", type=Path, default=DEFAULT_TLE_FILE)
    parser.add_argument("--start", default=DEFAULT_START_UTC.isoformat(), help="UTC start time in ISO-8601 format.")
    parser.add_argument("--seconds", type=float, default=90.0)
    parser.add_argument("--step-seconds", type=float, default=1.0)
    parser.add_argument(
        "--ephemeris-time-shift-seconds",
        type=float,
        default=TimeSeriesConfig.ephemeris_time_shift_seconds,
        help="Ephemeris time applied to the spacecraft state before plotting; positive values move the simulated pass earlier on the plotted UTC axis.",
    )
    parser.add_argument("--min-elevation-deg", type=float, default=0.0)
    parser.add_argument("--range-offset-km", type=float, default=0.0)
    parser.add_argument("--doppler-offset-hz", type=float, default=0.0)
    parser.add_argument("--reference-snr-db", type=float, default=SNRModel.reference_snr_db)
    parser.add_argument("--reference-range-km", type=float, default=SNRModel.reference_range_km)
    parser.add_argument("--reference-elevation-deg", type=float, default=SNRModel.reference_elevation_deg)
    parser.add_argument("--range-power-exponent", type=float, default=SNRModel.range_power_exponent)
    parser.add_argument("--beam-power-exponent", type=float, default=SNRModel.beam_power_exponent)
    parser.add_argument("--ripple-amplitude-db", type=float, default=SNRModel.ripple_amplitude_db)
    parser.add_argument("--ripple-period-km", type=float, default=SNRModel.ripple_period_km)
    parser.add_argument("--ripple-phase-rad", type=float, default=SNRModel.ripple_phase_rad)
    parser.add_argument("--floor-snr-db", type=float, default=SNRModel.floor_snr_db)
    parser.add_argument("--f107", type=float, default=None)
    parser.add_argument("--f107a", type=float, default=None)
    parser.add_argument("--measurement-csv", type=Path, default=DEFAULT_MEASUREMENT_FILE)
    parser.add_argument("--fit-geometry", action="store_true")
    parser.add_argument("--fit-snr", action="store_true")
    parser.add_argument("--fit-start", default=None, help="UTC fit start time in ISO-8601 format. Defaults to --start.")
    parser.add_argument("--fit-end", default=None, help="UTC fit end time in ISO-8601 format. Defaults to the end of the measurement file.")
    parser.add_argument("--hooke-amplitude-fraction", type=float, default=HookeWaveModel.amplitude_fraction)
    parser.add_argument("--hooke-horizontal-wavelength-km", type=float, default=HookeWaveModel.horizontal_wavelength_km)
    parser.add_argument("--hooke-bearing-deg", type=float, default=HookeWaveModel.bearing_deg)
    parser.add_argument("--hooke-bearing-offset-deg", type=float, default=HookeWaveModel.bearing_offset_deg)
    parser.add_argument("--hooke-period-seconds", type=float, default=HookeWaveModel.period_seconds)
    parser.add_argument("--hooke-phase-rad", type=float, default=HookeWaveModel.phase_rad)
    parser.add_argument("--hooke-vertical-center-km", type=float, default=HookeWaveModel.vertical_center_km)
    parser.add_argument("--hooke-vertical-sigma-km", type=float, default=HookeWaveModel.vertical_sigma_km)
    parser.add_argument("--hooke-snr-coupling-db", type=float, default=HookeWaveModel.snr_coupling_db)
    parser.add_argument("--raytube-delta-elevation-deg", type=float, default=RayTubeModel.delta_elevation_deg)
    parser.add_argument("--raytube-delta-bearing-deg", type=float, default=RayTubeModel.delta_bearing_deg)
    parser.add_argument("--real-grid-lat-step-deg", type=float, default=TimeSeriesConfig.real_grid_lat_step_deg)
    parser.add_argument("--real-grid-lon-step-deg", type=float, default=TimeSeriesConfig.real_grid_lon_step_deg)
    parser.add_argument("--real-grid-lat-margin-deg", type=float, default=TimeSeriesConfig.real_grid_lat_margin_deg)
    parser.add_argument("--real-grid-lon-margin-deg", type=float, default=TimeSeriesConfig.real_grid_lon_margin_deg)
    parser.add_argument("--real-grid-alt-min-km", type=float, default=TimeSeriesConfig.real_grid_alt_min_km)
    parser.add_argument("--real-grid-alt-max-km", type=float, default=TimeSeriesConfig.real_grid_alt_max_km)
    parser.add_argument("--real-grid-alt-step-km", type=float, default=TimeSeriesConfig.real_grid_alt_step_km)
    parser.add_argument("--real-d-region-model", choices=("fpt2018", "iri2020", "none"), default=TimeSeriesConfig.real_d_region_model)
    parser.add_argument("--csv-out", type=Path, default=Path("collm_skiymet_iss_prediction.csv"))
    parser.add_argument("--png-out", type=Path, default=Path("collm_skiymet_iss_prediction.png"))
    parser.add_argument("--fof2-map-out", type=Path, default=Path("collm_skiymet_iss_fof2_map.png"))
    parser.add_argument("--fof2-map-mode", choices=("total", "delta", "fractional"), default="total")
    parser.add_argument("--swath-out", type=Path, default=Path("collm_skiymet_iss_cpa_swath.png"))
    parser.add_argument("--fof2-lat-step-deg", type=float, default=0.25)
    parser.add_argument("--fof2-lon-step-deg", type=float, default=0.25)
    parser.add_argument("--fof2-margin-deg", type=float, default=6.0)
    parser.add_argument("--fof2-alt-min-km", type=float, default=150.0)
    parser.add_argument("--fof2-alt-max-km", type=float, default=500.0)
    parser.add_argument("--fof2-alt-step-km", type=float, default=10.0)
    parser.add_argument("--swath-waypoints", type=int, default=256)
    parser.add_argument("--swath-alt-min-km", type=float, default=0.0)
    parser.add_argument("--swath-alt-max-km", type=float, default=500.0)
    parser.add_argument("--swath-alt-step-km", type=float, default=5.0)
    parser.add_argument("--noise-seed", type=int, default=0)
    return parser


def _config_from_args(args: argparse.Namespace) -> TimeSeriesConfig:
    snr_model = SNRModel(
        reference_snr_db=float(args.reference_snr_db),
        reference_range_km=float(args.reference_range_km),
        reference_elevation_deg=float(args.reference_elevation_deg),
        range_power_exponent=float(args.range_power_exponent),
        beam_power_exponent=float(args.beam_power_exponent),
        ripple_amplitude_db=float(args.ripple_amplitude_db),
        ripple_period_km=float(args.ripple_period_km),
        ripple_phase_rad=float(args.ripple_phase_rad),
        floor_snr_db=float(args.floor_snr_db),
    )
    hooke_wave = HookeWaveModel(
        amplitude_fraction=float(args.hooke_amplitude_fraction),
        horizontal_wavelength_km=float(args.hooke_horizontal_wavelength_km),
        bearing_deg=None if args.hooke_bearing_deg is None else float(args.hooke_bearing_deg),
        bearing_offset_deg=float(args.hooke_bearing_offset_deg),
        period_seconds=float(args.hooke_period_seconds),
        phase_rad=float(args.hooke_phase_rad),
        vertical_center_km=float(args.hooke_vertical_center_km),
        vertical_sigma_km=float(args.hooke_vertical_sigma_km),
        snr_coupling_db=float(args.hooke_snr_coupling_db),
    )
    ray_tube_model = RayTubeModel(
        delta_elevation_deg=float(args.raytube_delta_elevation_deg),
        delta_bearing_deg=float(args.raytube_delta_bearing_deg),
    )
    return TimeSeriesConfig(
        tle_file=args.tle_file.expanduser(),
        start_utc=_parse_time(args.start),
        seconds=float(args.seconds),
        step_seconds=float(args.step_seconds),
        ephemeris_time_shift_seconds=float(args.ephemeris_time_shift_seconds),
        min_elevation_deg=float(args.min_elevation_deg),
        range_offset_km=float(args.range_offset_km),
        doppler_offset_hz=float(args.doppler_offset_hz),
        f107=None if args.f107 is None else float(args.f107),
        f107a=None if args.f107a is None else float(args.f107a),
        snr_model=snr_model,
        hooke_wave=hooke_wave,
        ray_tube_model=ray_tube_model,
        real_grid_lat_step_deg=float(args.real_grid_lat_step_deg),
        real_grid_lon_step_deg=float(args.real_grid_lon_step_deg),
        real_grid_lat_margin_deg=float(args.real_grid_lat_margin_deg),
        real_grid_lon_margin_deg=float(args.real_grid_lon_margin_deg),
        real_grid_alt_min_km=float(args.real_grid_alt_min_km),
        real_grid_alt_max_km=float(args.real_grid_alt_max_km),
        real_grid_alt_step_km=float(args.real_grid_alt_step_km),
        real_d_region_model=str(args.real_d_region_model),
        noise_seed=int(args.noise_seed),
    )


def main() -> None:
    args = build_parser().parse_args()
    config = _config_from_args(args)
    plot_measurements: tuple[MeasurementPoint, ...] = ()
    measurement_csv = None if args.measurement_csv is None else args.measurement_csv.expanduser()
    if measurement_csv is not None and measurement_csv.exists():
        plot_start = config.start_utc
        plot_end = config.start_utc + dt.timedelta(seconds=config.seconds)
        plot_measurements = load_measurements(measurement_csv, start_utc=plot_start, end_utc=plot_end)
    if args.fit_geometry:
        if measurement_csv is None or not measurement_csv.exists():
            raise SystemExit("measurement CSV is required for --fit-geometry")
        fit_start = config.start_utc if args.fit_start is None else _parse_time(args.fit_start)
        fit_end = None if args.fit_end is None else _parse_time(args.fit_end)
        fit_measurements = load_measurements(measurement_csv, start_utc=fit_start, end_utc=fit_end)
        fit = _fit_geometry_to_measurements(
            CollmSkiymetRadar(),
            config.tle_file,
            fit_measurements,
            initial_ephemeris_time_shift_seconds=config.ephemeris_time_shift_seconds,
            initial_range_offset_km=config.range_offset_km,
            initial_doppler_offset_hz=config.doppler_offset_hz,
        )
        config = replace(
            config,
            ephemeris_time_shift_seconds=fit.ephemeris_time_shift_seconds,
            range_offset_km=fit.range_offset_km,
            doppler_offset_hz=fit.doppler_offset_hz,
        )
        print(
            "Geometry fit: "
            f"time_shift={fit.ephemeris_time_shift_seconds:.3f} s, "
            f"range_offset={fit.range_offset_km:.3f} km, "
            f"doppler_offset={fit.doppler_offset_hz:.3f} Hz, "
            f"rms_range={fit.rms_range_km:.3f} km, "
            f"rms_doppler={fit.rms_doppler_hz:.3f} Hz"
        )
    try:
        series = predict_collm_skiymet_time_series(config)
    except PyLapImportError as exc:
        raise SystemExit(str(exc)) from exc
    if args.fit_snr:
        if measurement_csv is None or not measurement_csv.exists():
            raise SystemExit("measurement CSV is required for --fit-snr")
        fit_start = config.start_utc if args.fit_start is None else _parse_time(args.fit_start)
        fit_end = None if args.fit_end is None else _parse_time(args.fit_end)
        snr_measurements = load_measurements(measurement_csv, start_utc=fit_start, end_utc=fit_end)
        snr_fit = _fit_snr_to_measurements(series, snr_measurements)
        series = _apply_snr_fit(series, snr_fit)
        print(
            "SNR fit: "
            f"offset={snr_fit.offset_db:.3f} dB, "
            f"scale={snr_fit.scale:.3f}, "
            f"rms={snr_fit.rms_db:.3f} dB"
        )
    fof2_snapshot = build_fof2_map_snapshot(
        series,
        lat_step_deg=float(args.fof2_lat_step_deg),
        lon_step_deg=float(args.fof2_lon_step_deg),
        margin_deg=float(args.fof2_margin_deg),
        alt_min_km=float(args.fof2_alt_min_km),
        alt_max_km=float(args.fof2_alt_max_km),
        alt_step_km=float(args.fof2_alt_step_km),
    )
    swath_snapshot = build_cpa_swath_snapshot(
        series,
        waypoint_count=int(args.swath_waypoints),
        alt_min_km=float(args.swath_alt_min_km),
        alt_max_km=float(args.swath_alt_max_km),
        alt_step_km=float(args.swath_alt_step_km),
    )
    write_prediction_csv(args.csv_out.expanduser(), series)
    plot_prediction(args.png_out.expanduser(), series, measurements=plot_measurements)
    plot_fof2_map(
        args.fof2_map_out.expanduser(),
        fof2_snapshot,
        series.radar,
        mode=str(args.fof2_map_mode),
    )
    plot_cpa_swath(args.swath_out.expanduser(), swath_snapshot, series.radar)
    print(f"Wrote {args.csv_out.expanduser()}")
    print(f"Wrote {args.png_out.expanduser()}")
    print(f"Wrote {args.fof2_map_out.expanduser()}")
    print(f"Wrote {args.swath_out.expanduser()}")


if __name__ == "__main__":
    main()
