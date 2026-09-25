from __future__ import annotations

import datetime as dt
import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize

from .geometry import (
    GeoPoint,
    RayDistance,
    coerce_longitude_for_grid,
    initial_bearing_deg,
    llh_to_ecef,
    ray_point_distance,
    relaz_to_ecef_unit,
    vector_angle_deg,
    wrap_longitude,
)
from .grid import IonosphereGrid, build_pyiri_grid


@dataclass
class RayTrace:
    summary: dict
    path: dict
    state: dict
    home: bool = False
    error_m: float = math.inf
    group_range_to_rx_km: float | None = None
    geometric_dist_to_rx_km: float | None = None
    total_absorption_db: float | None = None
    perigee_km: float | None = None
    metadata: dict = field(default_factory=dict)

    @property
    def launch_elevation_deg(self) -> float:
        return float(self.path["initial_elev"])

    @property
    def launch_bearing_deg(self) -> float:
        return float(self.path["initial_bearing"])

    @property
    def frequency_mhz(self) -> float:
        return float(self.path["frequency"])

    @property
    def origin(self) -> GeoPoint:
        return GeoPoint(
            float(np.asarray(self.path["lat"])[0]),
            float(np.asarray(self.path["lon"])[0]),
            float(np.asarray(self.path["height"])[0]),
        )

    def to_dict(self) -> dict:
        return {
            "frequency_mhz": self.frequency_mhz,
            "home": self.home,
            "error_m": self.error_m,
            "launch_elevation_deg": self.launch_elevation_deg,
            "launch_bearing_deg": self.launch_bearing_deg,
            "group_range_to_rx_km": self.group_range_to_rx_km,
            "geometric_dist_to_rx_km": self.geometric_dist_to_rx_km,
            "total_absorption_db": self.total_absorption_db,
            "perigee_km": self.perigee_km,
        }


class PyLapImportError(RuntimeError):
    pass


def _ensure_pharlap_runtime_env() -> None:
    pharlap_home = Path(os.environ.get("PHARLAP_HOME", "/Users/chartat1/pharlap")).expanduser()
    if pharlap_home.is_dir():
        os.environ.setdefault("PHARLAP_HOME", str(pharlap_home))
        ref_dir = pharlap_home / "dat"
        if ref_dir.is_dir():
            os.environ.setdefault("DIR_MODELS_REF_DAT", str(ref_dir))


class PyLapRaytraceBackend:
    def __init__(self, raytrace_3d_func=None):
        if raytrace_3d_func is None:
            _ensure_pharlap_runtime_env()
            try:
                from pylap.raytrace_3d import raytrace_3d
            except ImportError as exc:
                raise PyLapImportError(
                    "pylap is not importable in this environment. Install or build pylap "
                    "before attempting a full raytrace."
                ) from exc
            raytrace_3d_func = raytrace_3d
        self._raytrace_3d = raytrace_3d_func

    def trace(self, origin: GeoPoint, elevations_deg: Sequence[float], bearings_deg: Sequence[float],
              freqs_mhz: Sequence[float], ox_mode: int, nhops: int, tol: Sequence[float], *,
              grid: IonosphereGrid | None = None, state_vector: dict[str, np.ndarray] | None = None) -> list[RayTrace]:
        elevs = np.asarray(elevations_deg, dtype=float)
        bearings = np.asarray(bearings_deg, dtype=float)
        freqs = np.asarray(freqs_mhz, dtype=float)
        args: list = [
            float(origin.lat_deg),
            float(wrap_longitude(origin.lon_deg)),
            float(origin.alt_km),
            elevs,
            bearings,
            freqs,
            int(ox_mode),
            int(nhops),
            [float(t) for t in tol],
        ]
        if grid is not None:
            args.extend([
                np.asarray(grid.iono_en_grid, dtype=float),
                np.asarray(grid.iono_en_grid_5, dtype=float),
                np.asarray(grid.collision_freq, dtype=float),
                [float(v) for v in grid.iono_grid_parms],
                np.asarray(grid.Bx, dtype=float),
                np.asarray(grid.By, dtype=float),
                np.asarray(grid.Bz, dtype=float),
                [float(v) for v in grid.geomag_grid_parms],
            ])
        if state_vector is not None:
            args.append({key: np.asarray(value, dtype=float) for key, value in state_vector.items()})

        ray_summaries, ray_paths, ray_states = self._raytrace_3d(*args)
        return [
            RayTrace(summary=summary, path=path, state=state)
            for summary, path, state in zip(ray_summaries, ray_paths, ray_states)
        ]


def _complex_to_valid_real(value: complex | float) -> float | None:
    if isinstance(value, complex):
        if abs(value.imag) > 1e-9:
            return None
        value = value.real
    value = float(value)
    if not math.isfinite(value) or value <= 0.0:
        return None
    return value


def appleton_hartree(theta_deg: float, ne_m3: float, b_tesla: float, wave_hz: float) -> tuple[float | None, float | None, float | None]:
    m_e = 9.10938356e-31
    charge_e = 1.60217662e-19
    epsilon_0 = 8.85418782e-12

    omega = 2.0 * math.pi * wave_hz
    omega_0 = math.sqrt((ne_m3 * charge_e ** 2) / (epsilon_0 * m_e))
    omega_h = abs(b_tesla) * abs(charge_e) / m_e

    x = (omega_0 ** 2) / (omega ** 2)
    y = omega_h / omega
    if abs(1.0 - x) < 1e-12:
        return None, None, None
    theta = math.radians(theta_deg)
    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)

    common_sqrt = math.sqrt(
        0.25 * y ** 4 * sin_theta ** 4 + y ** 2 * cos_theta ** 2 * (1.0 - x) ** 2
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        n_o = np.sqrt(
            1.0 - x / (
                1.0
                - 0.5 * y ** 2 * sin_theta ** 2 / (1.0 - x)
                + common_sqrt / (1.0 - x)
            )
        )
        n_x = np.sqrt(
            1.0 - x / (
                1.0
                - 0.5 * y ** 2 * sin_theta ** 2 / (1.0 - x)
                - common_sqrt / (1.0 - x)
            )
        )
        plasma_hz = omega_0 / (2.0 * math.pi)
        n_no_field = np.sqrt(1.0 - (plasma_hz ** 2) / (wave_hz ** 2))

    return (
        _complex_to_valid_real(n_o),
        _complex_to_valid_real(n_x),
        _complex_to_valid_real(n_no_field),
    )


class PointToPointRayTracer:
    def __init__(self, backend: PyLapRaytraceBackend | None = None):
        self.backend = backend

    def _backend(self) -> PyLapRaytraceBackend:
        if self.backend is None:
            self.backend = PyLapRaytraceBackend()
        return self.backend

    def prepare_transmitter_state(
        self, *, tx: GeoPoint, grid: IonosphereGrid,
    ) -> tuple[GeoPoint, float, float, float, float]:
        tx_local = self._nudge_if_on_grid(tx, grid)
        return (tx_local, *self._sample_tx_state(grid, tx_local))

    def prepare_ray_state_vector_batch(
        self,
        *,
        tx: GeoPoint,
        grid: IonosphereGrid,
        elevations_deg: Sequence[float],
        bearings_deg: Sequence[float],
        freqs_mhz: Sequence[float],
        ox_mode: int = 0,
        transmitter_state: tuple[GeoPoint, float, float, float, float] | None = None,
    ) -> tuple[GeoPoint, dict[str, np.ndarray], np.ndarray] | None:
        if transmitter_state is None:
            transmitter_state = self.prepare_transmitter_state(tx=tx, grid=grid)
        tx_local, ne_cm3, bx, by, bz = transmitter_state
        elevs = np.asarray(elevations_deg, dtype=float)
        bears = np.asarray(bearings_deg, dtype=float)
        freqs = np.asarray(freqs_mhz, dtype=float)
        state_vector = self._build_state_vector_batch(tx_local, elevs, bears, freqs, ox_mode, ne_cm3, bx, by, bz)
        if state_vector is None:
            return None
        valid_mask = np.asarray(state_vector.pop("_valid_mask"), dtype=bool)
        return tx_local, state_vector, valid_mask

    def trace_state_vector_batch(
        self,
        *,
        tx: GeoPoint,
        elevations_deg: Sequence[float],
        bearings_deg: Sequence[float],
        freqs_mhz: Sequence[float],
        grid: IonosphereGrid,
        state_vector: dict[str, np.ndarray],
        ox_mode: int = 0,
        nhops: int = 2,
        tol: Sequence[float] = (1e-8, 0.005, 5.0),
    ) -> list[RayTrace]:
        return self._backend().trace(
            tx,
            elevations_deg,
            bearings_deg,
            freqs_mhz,
            ox_mode,
            nhops,
            tol,
            grid=grid,
            state_vector=state_vector,
        )

    def trace_frequencies(self, *, when: dt.datetime, tx: GeoPoint, rx: GeoPoint,
                          frequencies_mhz: Sequence[float], f107: float | None = None,
                          grid: IonosphereGrid | None = None,
                          ox_mode: int = 0,
                          nhops: int = 2,
                          tol: Sequence[float] = (1e-8, 0.005, 5.0),
                          homing_tolerance_m: float = 100.0,
                          alt_min_km: float = 60.0,
                          alt_max_km: float = 500.0,
                          alt_step_km: float = 2.0,
                          lat_step_deg: float = 1.0,
                          lon_step_deg: float = 1.0,
                          lat_margin_deg: float = 5.0,
                          lon_margin_deg: float = 5.0,
                          d_region_model: str = "fpt2018",
                          ap_daily: float | None = None,
                          blend_bottom_km: float = 120.0,
                          blend_top_km: float = 140.0,
                          msis_version: float | str = 2.1,
                          refresh_indices: bool = False) -> list[RayTrace]:
        if grid is None:
            grid = build_pyiri_grid(
                when,
                tx,
                rx,
                f107=f107,
                alt_min_km=alt_min_km,
                alt_max_km=alt_max_km,
                alt_step_km=alt_step_km,
                lat_step_deg=lat_step_deg,
                lon_step_deg=lon_step_deg,
                lat_margin_deg=lat_margin_deg,
                lon_margin_deg=lon_margin_deg,
                d_region_model=d_region_model,
                ap_daily=ap_daily,
                blend_bottom_km=blend_bottom_km,
                blend_top_km=blend_top_km,
                msis_version=msis_version,
                refresh_indices=refresh_indices,
            )
        results = []
        for freq_mhz in frequencies_mhz:
            results.append(
                self.trace_link(
                    tx=tx,
                    rx=rx,
                    frequency_mhz=float(freq_mhz),
                    grid=grid,
                    ox_mode=ox_mode,
                    nhops=nhops,
                    tol=tol,
                    homing_tolerance_m=homing_tolerance_m,
                )
            )
        return results

    def trace_link(self, *, tx: GeoPoint, rx: GeoPoint, frequency_mhz: float, grid: IonosphereGrid,
                   ox_mode: int = 0, nhops: int = 2, tol: Sequence[float] = (1e-8, 0.005, 5.0),
                   homing_tolerance_m: float = 100.0,
                   start_angles_deg: tuple[float, float] | None = None) -> RayTrace:
        tx_local = self._nudge_if_on_grid(tx, grid)
        ne_cm3, bx, by, bz = self._sample_tx_state(grid, tx_local)
        best_start = self._initial_guess(
            tx=tx_local,
            rx=rx,
            frequency_mhz=frequency_mhz,
            grid=grid,
            ox_mode=ox_mode,
            nhops=nhops,
            tol=tol,
            ne_cm3=ne_cm3,
            bx=bx,
            by=by,
            bz=bz,
            start_angles_deg=start_angles_deg,
        )
        if best_start is None:
            return RayTrace(summary={}, path={"initial_elev": math.nan, "initial_bearing": math.nan, "frequency": frequency_mhz}, state={})

        return self._home_ray(
            start_ray=best_start,
            tx=tx_local,
            rx=rx,
            frequency_mhz=frequency_mhz,
            grid=grid,
            ox_mode=ox_mode,
            nhops=nhops,
            tol=tol,
            homing_tolerance_m=homing_tolerance_m,
            ne_cm3=ne_cm3,
            bx=bx,
            by=by,
            bz=bz,
        )

    def trace_launch(self, *, tx: GeoPoint, rx: GeoPoint, frequency_mhz: float, grid: IonosphereGrid,
                     elevation_deg: float, bearing_deg: float, ox_mode: int = 0, nhops: int = 2,
                     tol: Sequence[float] = (1e-8, 0.005, 5.0)) -> RayTrace:
        tx_local = self._nudge_if_on_grid(tx, grid)
        ne_cm3, bx, by, bz = self._sample_tx_state(grid, tx_local)
        return self._trace_candidate(
            tx=tx_local,
            rx=rx,
            grid=grid,
            elevation_deg=elevation_deg,
            bearing_deg=bearing_deg,
            frequency_mhz=frequency_mhz,
            ox_mode=ox_mode,
            nhops=nhops,
            tol=tol,
            ne_cm3=ne_cm3,
            bx=bx,
            by=by,
            bz=bz,
        )

    def _initial_guess(self, *, tx: GeoPoint, rx: GeoPoint, frequency_mhz: float, grid: IonosphereGrid,
                       ox_mode: int, nhops: int, tol: Sequence[float], ne_cm3: float, bx: float, by: float, bz: float,
                       start_angles_deg: tuple[float, float] | None) -> RayTrace | None:
        if start_angles_deg is not None:
            seeded = self._trace_candidate(
                tx=tx,
                rx=rx,
                grid=grid,
                elevation_deg=float(start_angles_deg[0]),
                bearing_deg=float(start_angles_deg[1]),
                frequency_mhz=frequency_mhz,
                ox_mode=ox_mode,
                nhops=nhops,
                tol=tol,
                ne_cm3=ne_cm3,
                bx=bx,
                by=by,
                bz=bz,
            )
            if math.isfinite(seeded.error_m):
                return seeded

        coarse_rays = self._global_search(
            tx=tx,
            rx=rx,
            frequency_mhz=frequency_mhz,
            grid=grid,
            ox_mode=ox_mode,
            nhops=nhops,
            tol=tol,
            ne_cm3=ne_cm3,
            bx=bx,
            by=by,
            bz=bz,
        )
        if not coarse_rays:
            return None
        return min(coarse_rays, key=lambda ray: ray.error_m)

    def _global_search(self, *, tx: GeoPoint, rx: GeoPoint, frequency_mhz: float, grid: IonosphereGrid,
                       ox_mode: int, nhops: int, tol: Sequence[float],
                       ne_cm3: float, bx: float, by: float, bz: float) -> list[RayTrace]:
        bearing_center = initial_bearing_deg(tx, rx)
        azimuths = np.arange(bearing_center - 3.0, bearing_center + 3.001, 1.0)
        elevations = np.arange(10.0, 90.001, 2.0)
        elevation_grid, azimuth_grid = np.meshgrid(elevations, azimuths, indexing="xy")
        elevs = elevation_grid.ravel()
        bears = azimuth_grid.ravel()
        freqs = np.full(elevs.shape, float(frequency_mhz), dtype=float)

        ray_state = self._build_state_vector_batch(tx, elevs, bears, freqs, ox_mode, ne_cm3, bx, by, bz)
        if ray_state is None:
            return []

        valid_mask = np.asarray(ray_state.pop("_valid_mask"), dtype=bool)
        rays = self._backend().trace(
            tx,
            elevs[valid_mask],
            bears[valid_mask],
            freqs[valid_mask],
            ox_mode,
            nhops,
            tol,
            grid=grid,
            state_vector=ray_state,
        )
        valid_rays: list[RayTrace] = []
        for ray in rays:
            distance = ray_point_distance(ray.path, rx)
            if math.isfinite(distance.distance_m):
                ray.error_m = distance.distance_m
                valid_rays.append(ray)
        return valid_rays

    def _home_ray(self, *, start_ray: RayTrace, tx: GeoPoint, rx: GeoPoint, frequency_mhz: float,
                  grid: IonosphereGrid, ox_mode: int, nhops: int, tol: Sequence[float], homing_tolerance_m: float,
                  ne_cm3: float, bx: float, by: float, bz: float) -> RayTrace:
        def objective(params: np.ndarray) -> float:
            result = self._trace_candidate(
                tx=tx,
                rx=rx,
                grid=grid,
                elevation_deg=float(params[0]),
                bearing_deg=float(params[1]),
                frequency_mhz=frequency_mhz,
                ox_mode=ox_mode,
                nhops=nhops,
                tol=tol,
                ne_cm3=ne_cm3,
                bx=bx,
                by=by,
                bz=bz,
            )
            return result.error_m

        x0 = np.array([start_ray.launch_elevation_deg, start_ray.launch_bearing_deg], dtype=float)
        result = minimize(
            objective,
            x0,
            method="Nelder-Mead",
            options={"fatol": homing_tolerance_m, "xatol": 0.005, "maxfev": 100, "disp": False},
        )
        homed = self._trace_candidate(
            tx=tx,
            rx=rx,
            grid=grid,
            elevation_deg=float(result.x[0]),
            bearing_deg=float(result.x[1]),
            frequency_mhz=frequency_mhz,
            ox_mode=ox_mode,
            nhops=nhops,
            tol=tol,
            ne_cm3=ne_cm3,
            bx=bx,
            by=by,
            bz=bz,
        )
        if homed.error_m < homing_tolerance_m:
            homed.home = True
            homed.perigee_km = float(np.nanmin(np.asarray(homed.path["height"], dtype=float)))
        return homed

    def _trace_candidate(self, *, tx: GeoPoint, rx: GeoPoint, grid: IonosphereGrid,
                         elevation_deg: float, bearing_deg: float, frequency_mhz: float,
                         ox_mode: int, nhops: int, tol: Sequence[float],
                         ne_cm3: float, bx: float, by: float, bz: float) -> RayTrace:
        if not (-90.0 <= elevation_deg <= 90.0):
            return RayTrace(summary={}, path={"initial_elev": elevation_deg, "initial_bearing": bearing_deg, "frequency": frequency_mhz}, state={})
        state_vector = self._build_state_vector_batch(
            tx,
            np.array([elevation_deg], dtype=float),
            np.array([bearing_deg], dtype=float),
            np.array([frequency_mhz], dtype=float),
            ox_mode,
            ne_cm3,
            bx,
            by,
            bz,
        )
        if state_vector is None:
            return RayTrace(summary={}, path={"initial_elev": elevation_deg, "initial_bearing": bearing_deg, "frequency": frequency_mhz}, state={})
        state_vector.pop("_valid_mask")
        try:
            ray = self._backend().trace(
                tx,
                [elevation_deg],
                [bearing_deg],
                [frequency_mhz],
                ox_mode,
                nhops,
                tol,
                grid=grid,
                state_vector=state_vector,
            )[0]
        except Exception:
            return RayTrace(summary={}, path={"initial_elev": elevation_deg, "initial_bearing": bearing_deg, "frequency": frequency_mhz}, state={})

        distance = ray_point_distance(ray.path, rx)
        ray.error_m = distance.distance_m
        if math.isfinite(distance.distance_m):
            ray.group_range_to_rx_km = distance.group_path_km
            ray.geometric_dist_to_rx_km = distance.geometric_path_km
            ray.total_absorption_db = distance.total_absorption_db
        return ray

    def _sample_tx_state(self, grid: IonosphereGrid, tx: GeoPoint) -> tuple[float, float, float, float]:
        if tx.alt_km < float(np.min(grid.altitudes_km)):
            return 0.0, 0.0, 0.0, 0.0

        lon = coerce_longitude_for_grid(tx.lon_deg, grid.longitudes_deg)
        point = np.array([[tx.lat_deg, lon, tx.alt_km]], dtype=float)

        density_interp = RegularGridInterpolator(
            (grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km),
            grid.iono_en_grid,
            bounds_error=False,
            fill_value=np.nan,
        )
        bx_interp = RegularGridInterpolator(
            (grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km),
            grid.Bx,
            bounds_error=False,
            fill_value=np.nan,
        )
        by_interp = RegularGridInterpolator(
            (grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km),
            grid.By,
            bounds_error=False,
            fill_value=np.nan,
        )
        bz_interp = RegularGridInterpolator(
            (grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km),
            grid.Bz,
            bounds_error=False,
            fill_value=np.nan,
        )

        ne_cm3 = float(density_interp(point)[0])
        bx = float(bx_interp(point)[0])
        by = float(by_interp(point)[0])
        bz = float(bz_interp(point)[0])
        if any(math.isnan(value) for value in (ne_cm3, bx, by, bz)):
            raise ValueError("failed to interpolate local transmitter state from the PyIRI grid")
        return ne_cm3, bx, by, bz

    def _build_state_vector_batch(self, tx: GeoPoint, elevations_deg: np.ndarray, bearings_deg: np.ndarray,
                                  freqs_mhz: np.ndarray, ox_mode: int, ne_cm3: float,
                                  bx: float, by: float, bz: float) -> dict[str, np.ndarray] | None:
        fields = {
            "pos_x": [],
            "pos_y": [],
            "pos_z": [],
            "dir_x": [],
            "dir_y": [],
            "dir_z": [],
            "group_path": [],
            "geometrical_path": [],
            "phase_path": [],
            "absorption": [],
            "indep_var": [],
            "ODE_step_size": [],
        }
        valid_mask: list[bool] = []

        pos = llh_to_ecef(tx.lat_deg, tx.lon_deg, tx.alt_km * 1000.0)
        if pos.ndim > 1:
            pos = pos[0]

        for elevation_deg, bearing_deg, freq_mhz in zip(elevations_deg, bearings_deg, freqs_mhz):
            state = self._calc_state_vector(
                pos_xyz_m=pos,
                tx=tx,
                elevation_deg=float(elevation_deg),
                bearing_deg=float(bearing_deg),
                freq_mhz=float(freq_mhz),
                ox_mode=ox_mode,
                ne_cm3=ne_cm3,
                bx=bx,
                by=by,
                bz=bz,
            )
            if state is None:
                valid_mask.append(False)
                continue
            valid_mask.append(True)
            for key, value in state.items():
                fields[key].append(value)

        if not any(valid_mask):
            return None

        out = {key: np.asarray(value, dtype=float) for key, value in fields.items()}
        out["_valid_mask"] = np.asarray(valid_mask, dtype=bool)
        return out

    def _calc_state_vector(self, *, pos_xyz_m: np.ndarray, tx: GeoPoint, elevation_deg: float, bearing_deg: float,
                           freq_mhz: float, ox_mode: int, ne_cm3: float, bx: float, by: float, bz: float) -> dict[str, float] | None:
        plasma_freq_mhz = math.sqrt(max(0.0, 80.6 * ne_cm3 / 1e6))
        if plasma_freq_mhz > freq_mhz:
            return None

        direction = np.asarray(relaz_to_ecef_unit(elevation_deg, bearing_deg, tx.lat_deg, tx.lon_deg), dtype=float)
        if direction.ndim > 1:
            direction = direction[0]

        b_vec = np.array([bx, by, bz], dtype=float)
        theta_deg = vector_angle_deg(direction, b_vec)
        b_mag = float(np.linalg.norm(b_vec))
        n_o, n_x, n_no_field = appleton_hartree(theta_deg, ne_cm3 * 1e6, b_mag, freq_mhz * 1e6)
        if ox_mode == 1:
            refractive_index = n_o
        elif ox_mode == -1:
            refractive_index = n_x
        else:
            refractive_index = n_no_field

        if refractive_index is None or refractive_index > 1.0:
            return None

        direction *= refractive_index
        return {
            "pos_x": float(pos_xyz_m[0]),
            "pos_y": float(pos_xyz_m[1]),
            "pos_z": float(pos_xyz_m[2]),
            "dir_x": float(direction[0]),
            "dir_y": float(direction[1]),
            "dir_z": float(direction[2]),
            "group_path": 0.0,
            "geometrical_path": 0.0,
            "phase_path": 0.0,
            "absorption": 0.0,
            "indep_var": 0.0,
            "ODE_step_size": 1000.0,
        }

    def _nudge_if_on_grid(self, tx: GeoPoint, grid: IonosphereGrid) -> GeoPoint:
        lat = tx.lat_deg + (1e-6 if np.any(np.isclose(grid.latitudes_deg, tx.lat_deg)) else 0.0)
        lon = tx.lon_deg + (1e-6 if np.any(np.isclose(grid.longitudes_deg, tx.lon_deg)) else 0.0)
        return GeoPoint(lat, lon, tx.alt_km)
