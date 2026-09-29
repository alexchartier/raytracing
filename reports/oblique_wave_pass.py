"""Simulate and fit a 600 km, two-satellite link across the saved wave pass.

Commands: trace, batch, build, select, evaluate. The fit uses only oblique O/X
ionograms until selection is frozen; truth density is opened by evaluate.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import resource
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator, UnivariateSpline
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "reports") not in sys.path:
    sys.path.insert(0, str(ROOT / "reports"))

from python_raytrace.absorption import effective_collision_frequency  # noqa: E402
from python_raytrace.geometry import GeoPoint, llh_to_ecef  # noqa: E402
from python_raytrace.grid import load_ionosphere_grid_netcdf, save_ionosphere_grid_netcdf  # noqa: E402
from python_raytrace.multisat_topside_inverse_demo import (  # noqa: E402
    TopsideInverseConfig, _fan_mesh, _topside_search_arrays,
    home_frequency_sweep_adaptive,
)
from python_raytrace.spacecraft_doppler import spacecraft_doppler  # noqa: E402
from ionogram_metrics import Ionogram, score  # noqa: E402
from build_lat_wave_doppler_peak_round3 import basis  # noqa: E402
from local_ray_lock import local_ray_lock, require_remote_for_large_grid  # noqa: E402

DATA = ROOT / "reports/data"
CASE_ROOT = Path(os.environ.get("SOUNDER_CASE_ROOT", str(DATA))).resolve()
RUN = CASE_ROOT / "oblique_wave_600km"
MANIFEST = CASE_ROOT / ("manifest.json" if CASE_ROOT != DATA else "lat_wave_pass_manifest.json")
TRUTH_GRID = CASE_ROOT / ("truth_grid.nc" if CASE_ROOT != DATA else "lat_wave_pass_forward_grid.nc")
TRUTH_DENSITY = CASE_ROOT / ("truth_density.npz" if CASE_ROOT != DATA else "lat_wave_pass_truth_density.npz")
START_GRID = Path(os.environ.get("SOUNDER_START_GRID", str(DATA / "lat_wave_doppler_peak_round3/wide_grid.nc"))).resolve()
START_DENSITY = Path(os.environ.get("SOUNDER_START_DENSITY", str(DATA / "lat_wave_doppler_peak_round3/wide_density.npz"))).resolve()
FIT_METADATA = CASE_ROOT / ("vertical_fit.json" if CASE_ROOT != DATA else "lat_wave_pass_initial_fit.json")
FIGURES = ROOT / "reports/figures"
FREQUENCIES = np.round(np.arange(2.0, 10.0001, 0.1), 10)
SEPARATION_KM = 600.0
SPACECRAFT_SPEED_MPS = 8000.0
RANGE_MIN_KM = 700.0


def profiles() -> list[dict]:
    return json.loads(MANIFEST.read_text())["profiles"]


def endpoints(index: int) -> tuple[GeoPoint, GeoPoint]:
    row = profiles()[index - 1]
    tx = GeoPoint(float(row["latitude_deg"]), float(row["longitude_deg"]),
                  float(row["altitude_km"]))
    def residual(latitude_offset: float) -> float:
        point = GeoPoint(tx.lat_deg + latitude_offset, tx.lon_deg, tx.alt_km)
        return separation_km(tx, point) - SEPARATION_KM

    latitude_offset = brentq(residual, 0.1, 10.0)
    rx = GeoPoint(tx.lat_deg + latitude_offset, tx.lon_deg, tx.alt_km)
    return tx, rx


def separation_km(tx: GeoPoint, rx: GeoPoint) -> float:
    a = np.asarray(llh_to_ecef(tx.lat_deg, tx.lon_deg, tx.alt_km * 1000.0))
    b = np.asarray(llh_to_ecef(rx.lat_deg, rx.lon_deg, rx.alt_km * 1000.0))
    return float(np.linalg.norm(a - b) / 1000.0)


def trace(grid_path: Path, index: int, output: Path, density_source: str) -> dict:
    started = time.perf_counter()
    require_remote_for_large_grid(grid_path)
    tx, rx = endpoints(index)
    if not math.isclose(separation_km(tx, rx), SEPARATION_KM, abs_tol=0.01):
        raise ValueError("Satellite spacing is not 600 km")
    grid = load_ionosphere_grid_netcdf(grid_path)
    if not (grid.latitudes_deg[0] < tx.lat_deg < rx.lat_deg < grid.latitudes_deg[-1]
            and grid.longitudes_deg[0] < tx.lon_deg < grid.longitudes_deg[-1]):
        raise ValueError("Both satellites must lie inside the forward grid")
    config = replace(TopsideInverseConfig(), oblique_elevation_count=15,
                     oblique_bearing_count=5, seed_max_candidates_per_frequency=32,
                     homed_max_returns_per_frequency=64, d_region_model="none")
    elevations, bearings = _fan_mesh(*_topside_search_arrays(tx, rx, config))
    records: list[tuple[float, ...]] = []
    dopplers: list[float] = []
    launch: list[tuple[float, float]] = []
    arrival: list[tuple[float, float]] = []
    perigees: list[float] = []
    counts = np.zeros((len(FREQUENCIES), 2), dtype=int)
    for mode_column, mode in enumerate((1, -1)):
        sweep = home_frequency_sweep_adaptive(
            tx=tx, rx=rx, grid=grid,
            fan_elevations_deg=elevations, fan_bearings_deg=bearings,
            frequencies_mhz=FREQUENCIES, ox_mode=mode, config=config,
            range_min_km=RANGE_MIN_KM, anchor_stride=5, anchor_block_size=10,
            anchor_elevation_stride=1, anchor_optimizer="Nelder-Mead",
        )
        for frequency_index, returns in enumerate(sweep):
            counts[frequency_index, mode_column] = len(returns)
            for found in returns:
                minimum_height = float(np.min(np.asarray(found.ray.path["height"], dtype=float)))
                if minimum_height >= min(tx.alt_km, rx.alt_km) - 100.0:
                    raise ValueError("Direct ray passed the reflected-return check")
                observable = spacecraft_doppler(
                    found.ray, tx, rx, float(FREQUENCIES[frequency_index]),
                    speed_mps=SPACECRAFT_SPEED_MPS, track_bearing_deg=0.0)
                records.append((float(frequency_index), float(mode),
                                found.group_range_km, found.miss_m,
                                found.absorption_db))
                dopplers.append(observable.doppler_hz)
                launch.append((observable.launch_elevation_deg,
                               observable.launch_bearing_deg))
                arrival.append((observable.arrival_elevation_deg,
                                observable.arrival_bearing_deg))
                perigees.append(minimum_height)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output, records=np.asarray(records, dtype=float).reshape(-1, 5),
        frequencies_mhz=FREQUENCIES, count_array=counts,
        spacecraft_doppler_hz=np.asarray(dopplers),
        launch_angles_deg=np.asarray(launch, dtype=float).reshape(-1, 2),
        arrival_angles_deg=np.asarray(arrival, dtype=float).reshape(-1, 2),
        ray_minimum_altitude_km=np.asarray(perigees),
        tx_lat_deg=np.array(tx.lat_deg), tx_lon_deg=np.array(tx.lon_deg),
        tx_alt_km=np.array(tx.alt_km), rx_lat_deg=np.array(rx.lat_deg),
        rx_lon_deg=np.array(rx.lon_deg), rx_alt_km=np.array(rx.alt_km),
        satellite_separation_km=np.array(separation_km(tx, rx)),
        spacecraft_speed_mps=np.array(SPACECRAFT_SPEED_MPS),
        doppler_model=np.array("local transmitter and receiver tangent velocities"),
        profile_index=np.array(index), density_source=np.array(density_source),
        method=np.array("oblique_adaptive"),
        vertical_fan_layout=np.array("along_track"),
        vertical_outer_ray_fraction=np.array(1.0),
        vertical_guard_seed_limit=np.array(-1),
        fan_launch_directions=np.array(len(elevations)),
        homing_tolerance_m=np.array(config.homing_tolerance_m),
        runtime_seconds=np.array(time.perf_counter() - started),
        peak_rss_mb=np.array(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0),
    )
    return {"profile": index, "accepted_returns": len(records),
            "runtime_seconds": time.perf_counter() - started,
            "nose_mhz_by_mode": [float(FREQUENCIES[np.max(np.flatnonzero(counts[:, col]))])
                                 if np.any(counts[:, col]) else None for col in (0, 1)]}


def batch(grid: Path, output_dir: Path, workers: int, indices: list[int],
          density_source: str) -> dict:
    if workers != 1:
        raise ValueError("local memory safety requires exactly one ray worker")
    if not indices or any(index < 1 or index > 20 for index in indices):
        raise ValueError("profile index must be in 1..20")
    output_dir.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
    started = time.perf_counter()

    def one(index: int) -> dict:
        destination = output_dir / f"ionogram_{index:02d}.npz"
        if destination.exists():
            with np.load(destination, allow_pickle=False) as saved:
                if (int(saved["profile_index"]) == index
                        and str(saved["method"]) == "oblique_adaptive"
                        and "doppler_model" in saved
                        and str(saved["doppler_model"]) == "local transmitter and receiver tangent velocities"
                        and len(saved["frequencies_mhz"]) == 81
                        and len(saved["spacecraft_doppler_hz"]) == len(saved["records"])):
                    return {"profile": index, "accepted_returns": len(saved["records"]),
                            "runtime_seconds": float(saved["runtime_seconds"]), "reused": True}
        command = [sys.executable, str(Path(__file__).resolve()), "trace",
                   "--grid", str(grid), "--index", str(index), "--output", str(destination),
                   "--density-source", density_source]
        result = subprocess.run(command, cwd=ROOT, env=environment,
                                capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(f"Profile {index} failed:\n{result.stdout}\n{result.stderr}")
        return {**json.loads(result.stdout), "reused": False}

    completed = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one, index): index for index in indices}
        for future in as_completed(futures):
            item = future.result()
            completed.append(item)
            print(f"{output_dir.name}: profile {item['profile']:02d}, "
                  f"{item['accepted_returns']} returns", flush=True)
    summary = {"grid": str(grid.resolve().relative_to(ROOT)),
               "output_dir": str(output_dir.resolve().relative_to(ROOT)),
               "workers": workers, "elapsed_seconds": time.perf_counter() - started,
               "profiles": sorted(completed, key=lambda row: row["profile"])}
    (output_dir / "batch.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def build_candidate(gain: float, vertical_sigma_km: float,
                    height_gain: float = 0.0, height_fit: str = "basis") -> dict:
    """Infer smooth peak and reflector-height corrections from oblique O/X rays."""
    if not 0.0 <= gain <= 1.5 or not 0.0 <= height_gain <= 2.0:
        raise ValueError("Invalid correction gain")
    if gain == 0.0 and height_gain == 0.0:
        raise ValueError("At least one correction must be nonzero")
    if not 20.0 <= vertical_sigma_km <= 250.0:
        raise ValueError("Invalid vertical width")
    if height_fit not in ("basis", "spline"):
        raise ValueError("height fit must be basis or spline")
    if not (RUN / "truth").is_dir() or not (RUN / "baseline").is_dir():
        raise FileNotFoundError("Trace truth and baseline oblique passes first")
    tx_latitudes = np.array([row["latitude_deg"] for row in profiles()], dtype=float)
    mid_latitudes = np.array([np.mean([endpoints(i)[0].lat_deg,
                                       endpoints(i)[1].lat_deg]) for i in range(1, 21)])
    center = float(json.loads(MANIFEST.read_text())["latitude_center_deg"])
    halfspan_km = float(np.max(abs(6371.0088 * np.deg2rad(tx_latitudes - center))))
    wavelength_km = float(json.loads(FIT_METADATA.read_text())["wave_wavelength_km"])
    design = basis(mid_latitudes, center, halfspan_km, wavelength_km)
    available = sorted(int(path.stem.split("_")[-1]) for path in (RUN / "truth").glob("ionogram_*.npz")
                       if (RUN / "baseline" / path.name).is_file())
    if len(available) < 4:
        raise ValueError("At least four paired oblique ionograms are required")
    noses = np.full((20, 2, 2), np.nan)
    low_frequency_range_residual_km = np.full(20, np.nan)
    for index in available:
        observed_ionogram = Ionogram.read(RUN / "truth" / f"ionogram_{index:02d}.npz")
        modeled_ionogram = Ionogram.read(RUN / "baseline" / f"ionogram_{index:02d}.npz")
        for source_column, source in enumerate(("truth", "baseline")):
            ionogram = observed_ionogram if source == "truth" else modeled_ionogram
            for mode_column, mode in enumerate((1, -1)):
                nose = ionogram.nose(mode)
                if nose is not None:
                    noses[index - 1, source_column, mode_column] = nose
        residuals = []
        for mode in (1, -1):
            truth_ridge = observed_ionogram.ridge(mode)
            model_ridge = modeled_ionogram.ridge(mode)
            for frequency_index in set(truth_ridge) & set(model_ridge):
                frequency = FREQUENCIES[frequency_index]
                if 2.5 <= frequency <= 4.5:
                    residuals.append(truth_ridge[frequency_index]
                                     - model_ridge[frequency_index])
        if residuals:
            low_frequency_range_residual_km[index - 1] = float(np.median(residuals))
    observed, modeled = noses[:, 0], noses[:, 1]
    valid = np.isfinite(observed) & np.isfinite(modeled) & (modeled < 9.9)
    proxy = np.clip(2.0 * (observed - modeled) / modeled, -.20, .20)
    rows = np.flatnonzero(np.any(valid, axis=1))
    if len(rows) < 4:
        raise ValueError("Too few paired oblique cutoffs for a wave fit")
    target = np.array([np.mean(proxy[i, valid[i]]) for i in rows])
    # Soft prior limits extrapolation at the ends of the pass.
    ridge = np.diag([.8, 1.2, 1.2, 1.2])
    coefficients = np.linalg.lstsq(
        np.vstack((design[rows], ridge)),
        np.r_[target, np.zeros(4)], rcond=None)[0]
    height_rows = np.flatnonzero(np.isfinite(low_frequency_range_residual_km))
    if len(height_rows) < 4:
        raise ValueError("Too few common low-frequency ranges for height fit")
    height_target_km = -low_frequency_range_residual_km[height_rows] / 1.7
    height_coefficients = np.linalg.lstsq(
        np.vstack((design[height_rows], ridge)),
        np.r_[height_target_km, np.zeros(4)], rcond=None)[0]
    grid = load_ionosphere_grid_netcdf(START_GRID)
    density = np.asarray(grid.iono_en_grid, dtype=float)
    grid_design = basis(np.asarray(grid.latitudes_deg), center,
                        halfspan_km, wavelength_km)
    fraction = gain * (grid_design @ coefficients)
    if height_fit == "spline":
        height_spline = UnivariateSpline(mid_latitudes[height_rows],
                                        height_target_km, s=20.0, k=3)
        clipped_latitudes = np.clip(np.asarray(grid.latitudes_deg),
                                    mid_latitudes[height_rows[0]],
                                    mid_latitudes[height_rows[-1]])
        height_shift_km = height_gain * height_spline(clipped_latitudes)
    else:
        height_shift_km = height_gain * (grid_design @ height_coefficients)
    if np.max(abs(fraction)) > .12:
        raise ValueError("Oblique correction exceeds 12% trust region")
    if np.max(abs(height_shift_km)) > 25.0:
        raise ValueError("Oblique height correction exceeds 25 km trust region")
    altitude = np.asarray(grid.altitudes_km, dtype=float)
    shifted = density.copy()
    for latitude_index, shift in enumerate(height_shift_km):
        for longitude_index in range(density.shape[1]):
            shifted[latitude_index, longitude_index] = np.interp(
                altitude - shift, altitude, density[latitude_index, longitude_index],
                left=density[latitude_index, longitude_index, 0],
                right=density[latitude_index, longitude_index, -1])
    peak_altitude = altitude[np.argmax(shifted, axis=2)]
    envelope = np.exp(-.5 * ((altitude[None, None, :]
                             - peak_altitude[:, :, None]) / vertical_sigma_km) ** 2)
    revised = shifted * (1.0 + fraction[:, None, None] * envelope)
    collision = grid.collision_freq
    if (grid.electron_temp_k is not None and grid.ion_temp_k is not None
            and grid.neutral_species_cm3 is not None):
        collision = effective_collision_frequency(
            grid.electron_temp_k, grid.ion_temp_k, revised * 1e6,
            grid.neutral_species_cm3)
    name = f"gain{gain:g}_sigma{vertical_sigma_km:g}"
    if height_gain:
        name += f"_height{height_gain:g}"
        if height_fit == "spline":
            name += "_spline"
    directory = RUN / name
    directory.mkdir(parents=True, exist_ok=True)
    grid_path = directory / "grid.nc"
    density_path = directory / "density.npz"
    save_ionosphere_grid_netcdf(
        grid_path, replace(grid, iono_en_grid=revised,
                           iono_en_grid_5=revised, collision_freq=collision))
    np.savez_compressed(density_path, latitudes_deg=grid.latitudes_deg,
                        longitudes_deg=grid.longitudes_deg,
                        altitudes_km=grid.altitudes_km,
                        electron_density_cm3=revised,
                        model=np.array("PyIRI vertical retrieval plus oblique cutoff fit"))
    result = {"selection_uses_truth_density": False,
              "observed_ionograms": str((RUN / "truth").relative_to(ROOT)),
              "modeled_ionograms": str((RUN / "baseline").relative_to(ROOT)),
              "starting_grid": str(START_GRID.relative_to(ROOT)),
              "satellite_separation_km": SEPARATION_KM,
              "basis": "constant, latitude trend, sine/cosine at prior fitted wavelength",
              "wavelength_km": wavelength_km,
              "coefficients": coefficients.tolist(),
              "height_coefficients_km": height_coefficients.tolist(),
              "low_frequency_range_residual_km": low_frequency_range_residual_km.tolist(),
              "height_gain": height_gain, "gain": gain,
              "height_fit": height_fit,
              "vertical_sigma_km": vertical_sigma_km,
              "height_shift_range_km": [float(height_shift_km.min()),
                                        float(height_shift_km.max())],
              "correction_range_percent": (100 * np.array([fraction.min(), fraction.max()])).tolist(),
              "paired_noses_mhz": noses.tolist(),
              "grid": str(grid_path.relative_to(ROOT)),
              "density": str(density_path.relative_to(ROOT))}
    (directory / "candidate.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def doppler_score(observed_path: Path, modeled_path: Path) -> float:
    """Use 1 Hz Doppler after matching each accepted return in mode/frequency/range."""
    with np.load(observed_path, allow_pickle=False) as source:
        if ("doppler_model" not in source or
                str(source["doppler_model"]) != "local transmitter and receiver tangent velocities"):
            raise ValueError(f"Stale Doppler model in {observed_path}")
        observed = np.asarray(source["records"], dtype=float)
        observed_doppler = np.rint(np.asarray(source["spacecraft_doppler_hz"], dtype=float))
    with np.load(modeled_path, allow_pickle=False) as source:
        if ("doppler_model" not in source or
                str(source["doppler_model"]) != "local transmitter and receiver tangent velocities"):
            raise ValueError(f"Stale Doppler model in {modeled_path}")
        modeled = np.asarray(source["records"], dtype=float)
        modeled_doppler = np.rint(np.asarray(source["spacecraft_doppler_hz"], dtype=float))
    errors = []
    for index, row in enumerate(observed):
        matches = np.flatnonzero((modeled[:, 0] == row[0]) & (modeled[:, 1] == row[1]))
        if not len(matches):
            errors.append(1.0)
            continue
        nearest = matches[np.argmin(abs(modeled[matches, 2] - row[2]))]
        if abs(modeled[nearest, 2] - row[2]) > 40.0:
            errors.append(1.0)
        else:
            errors.append(min(abs(observed_doppler[index] - modeled_doppler[nearest]) / 20.0, 1.0))
    return float(np.mean(errors)) if errors else 1.0


def select(names: list[str], indices: list[int], output: Path | None = None,
           doppler_weight: float = .2) -> dict:
    if not names or any(not (RUN / name / "ionograms").is_dir() for name in names):
        raise ValueError("Trace each named candidate before selection")
    if not indices or any(index < 1 or index > 20 for index in indices):
        raise ValueError("Profile indices must be in 1..20")
    indices = sorted(set(indices))
    output = RUN / "selection.json" if output is None else output
    if output == RUN / "selection.json" and indices != list(range(1, 21)):
        raise ValueError("Final selection requires all 20 profiles")
    rows = []
    for name in ["baseline", *names]:
        directory = RUN / ("baseline" if name == "baseline" else f"{name}/ionograms")
        ionogram_scores = []
        doppler_scores = []
        for index in indices:
            filename = f"ionogram_{index:02d}.npz"
            observed = RUN / "truth" / filename
            predicted = directory / filename
            ionogram_scores.append(score(Ionogram.read(observed), Ionogram.read(predicted)))
            if doppler_weight:
                doppler_scores.append(doppler_score(observed, predicted))
        mean_ionogram = float(np.mean([item["total"] for item in ionogram_scores]))
        mean_doppler = float(np.mean(doppler_scores)) if doppler_scores else None
        rows.append({"name": name, "ionograms": str(directory.relative_to(ROOT)),
                     "mean_ionogram_score": mean_ionogram,
                     "mean_1hz_doppler_score": mean_doppler,
                     "combined_score": ((1 - doppler_weight) * mean_ionogram
                                        + doppler_weight * mean_doppler if doppler_scores
                                        else mean_ionogram),
                     "ionogram_scores_by_profile": ionogram_scores,
                     "doppler_scores_by_profile": doppler_scores})
    chosen = min(rows, key=lambda row: row["combined_score"])
    result = {"selection_uses_truth_density": False,
              "selection_rule": ("minimum 0.8 O/X ionogram + 0.2 1-Hz Doppler score"
                                 if doppler_weight else "minimum O/X ionogram score"),
              "stage": "final" if indices == list(range(1, 21)) else "pilot",
              "profiles": indices, "chosen_candidate": chosen["name"],
              "candidates": rows}
    output.write_text(json.dumps(result, indent=2) + "\n")
    return result


def density_at_points(path: Path, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as source:
        altitudes = np.asarray(source["altitudes_km"], dtype=float)
        interpolator = RegularGridInterpolator(
            (source["latitudes_deg"], source["longitudes_deg"]),
            source["electron_density_cm3"], bounds_error=True)
        return altitudes, np.asarray(interpolator(points), dtype=float)


def evaluate() -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    selection = json.loads((RUN / "selection.json").read_text())
    if selection["selection_uses_truth_density"] or selection["profiles"] != list(range(1, 21)):
        raise ValueError("Expected a frozen 20-profile ionogram-only selection")
    name = selection["chosen_candidate"]
    fitted_dir = RUN / ("baseline" if name == "baseline" else f"{name}/ionograms")
    fitted_density = START_DENSITY if name == "baseline" else RUN / name / "density.npz"
    tx_points = np.array([[p["latitude_deg"], p["longitude_deg"]] for p in profiles()])
    midpoint_points = np.array([[(endpoints(i)[0].lat_deg + endpoints(i)[1].lat_deg) / 2,
                                 endpoints(i)[0].lon_deg] for i in range(1, 21)])
    metrics = {}
    sections = {}
    for label, points in (("link_midpoint", midpoint_points), ("transmitter", tx_points)):
        altitude, truth = density_at_points(TRUTH_DENSITY, points)
        start_altitude, starting = density_at_points(START_DENSITY, points)
        fit_altitude, fitted = density_at_points(fitted_density, points)
        if not (np.array_equal(altitude, start_altitude)
                and np.array_equal(altitude, fit_altitude)):
            raise ValueError("Truth and retrieved altitude axes differ")
        true_peak = np.max(truth, axis=1)
        starting_peak = np.max(starting, axis=1)
        fitted_peak = np.max(fitted, axis=1)
        mask = (altitude >= 150.0) & (altitude <= 600.0)
        errors = 100.0 * (fitted_peak / true_peak - 1.0)
        metrics[label] = {
            "baseline_peak_density_mae_percent": float(np.mean(abs(100 * (starting_peak / true_peak - 1)))),
            "retrieved_peak_density_mae_percent": float(np.mean(abs(errors))),
            "retrieved_peak_density_max_absolute_error_percent": float(np.max(abs(errors))),
            "baseline_profile_mean_nrmse_150_to_600_km": float(np.mean(
                np.sqrt(np.mean((starting[:, mask] - truth[:, mask]) ** 2, axis=1)) / true_peak)),
            "retrieved_profile_mean_nrmse_150_to_600_km": float(np.mean(
                np.sqrt(np.mean((fitted[:, mask] - truth[:, mask]) ** 2, axis=1)) / true_peak)),
            "retrieved_peak_error_percent_by_profile": errors.tolist(),
        }
        sections[label] = (altitude, truth, starting, fitted, true_peak)
    chosen = next(row for row in selection["candidates"] if row["name"] == name)
    baseline = next(row for row in selection["candidates"] if row["name"] == "baseline")
    midpoint_altitude, midpoint_truth, _, _, midpoint_true_peak = sections["link_midpoint"]
    midpoint_mask = (midpoint_altitude >= 150.0) & (midpoint_altitude <= 600.0)
    other_candidates = []
    for row in selection["candidates"]:
        other_name = row["name"]
        if other_name in ("baseline", name):
            continue
        _, alternate_density = density_at_points(RUN / other_name / "density.npz",
                                                  midpoint_points)
        alternate_peak = np.max(alternate_density, axis=1)
        other_candidates.append({
            "name": other_name,
            "combined_observable_score": row["combined_score"],
            "link_midpoint_peak_density_mae_percent": float(np.mean(abs(
                100.0 * (alternate_peak / midpoint_true_peak - 1.0)))),
            "link_midpoint_profile_mean_nrmse_150_to_600_km": float(np.mean(
                np.sqrt(np.mean((alternate_density[:, midpoint_mask]
                                 - midpoint_truth[:, midpoint_mask]) ** 2,
                                axis=1)) / midpoint_true_peak)),
        })
    observed_doppler = []
    observed_launch_elevation = []
    observed_arrival_elevation = []
    for index in range(1, 21):
        with np.load(RUN / "truth" / f"ionogram_{index:02d}.npz", allow_pickle=False) as source:
            observed_doppler.extend(np.asarray(source["spacecraft_doppler_hz"], dtype=float))
            observed_launch_elevation.extend(np.asarray(source["launch_angles_deg"], dtype=float)[:, 0])
            observed_arrival_elevation.extend(np.asarray(source["arrival_angles_deg"], dtype=float)[:, 0])
    result = {"selection_uses_truth_density": False,
              "evaluation_uses_truth_density_after_selection": True,
              "prior_vertical_retrieval_used_same_synthetic_truth_case": True,
              "truth_density_source": "Fortran IRI-2016 with imposed 900 km wave",
              "retrieval_start": str(START_DENSITY.relative_to(ROOT)),
              "selected_candidate": name,
              "satellite_spacing_km": SEPARATION_KM,
              "modeled_spacecraft_speed_mps": SPACECRAFT_SPEED_MPS,
              "homing_tolerance_m": 1000.0,
              "frequency_range_mhz": [2.0, 10.0], "frequency_step_mhz": 0.1,
              "truth_doppler_range_hz": [float(np.min(observed_doppler)),
                                         float(np.max(observed_doppler))],
              "truth_launch_elevation_range_deg": [float(np.min(observed_launch_elevation)),
                                                   float(np.max(observed_launch_elevation))],
              "truth_arrival_elevation_range_deg": [float(np.min(observed_arrival_elevation)),
                                                    float(np.max(observed_arrival_elevation))],
              "truth_accepted_returns": int(sum(len(Ionogram.read(
                  RUN / "truth" / f"ionogram_{i:02d}.npz").records) for i in range(1, 21))),
              "retrieved_accepted_returns": int(sum(len(Ionogram.read(
                  fitted_dir / f"ionogram_{i:02d}.npz").records) for i in range(1, 21))),
              "baseline_ionogram_score": baseline["mean_ionogram_score"],
              "retrieved_ionogram_score": chosen["mean_ionogram_score"],
              "baseline_doppler_score": baseline["mean_1hz_doppler_score"],
              "retrieved_doppler_score": chosen["mean_1hz_doppler_score"],
              "baseline_combined_score": baseline["combined_score"],
              "retrieved_combined_score": chosen["combined_score"],
              "other_candidates_after_selection": other_candidates,
              "density_metrics": metrics}
    (RUN / "evaluation.json").write_text(json.dumps(result, indent=2) + "\n")
    FIGURES.mkdir(parents=True, exist_ok=True)

    indices = (3, 10, 14, 18)
    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex=True, sharey=True,
                             constrained_layout=True)
    panel_maximum = 0.0
    for row, index in enumerate(indices):
        for col, (directory, label) in enumerate(((RUN / "truth", "Truth"),
                                                  (fitted_dir, "Retrieved"))):
            ax = axes[row, col]
            with np.load(directory / f"ionogram_{index:02d}.npz", allow_pickle=False) as source:
                records = np.asarray(source["records"], dtype=float)
                frequencies = np.asarray(source["frequencies_mhz"], dtype=float)
            if len(records):
                panel_maximum = max(panel_maximum, float(np.max(records[:, 2])))
            for mode, mode_label, color in ((1, "O", "#135b9a"),
                                            (-1, "X", "#b63836")):
                rays = records[records[:, 1] == mode]
                ax.scatter(frequencies[rays[:, 0].astype(int)], np.round(rays[:, 2]),
                           s=7, marker="s", linewidths=0, alpha=.8,
                           color=color, label=mode_label)
            ax.set(xlim=(2, 10),
                   title=f"{label}, link {index:02d} ({midpoint_points[index-1, 0]:.1f}° midpoint)")
            ax.grid(alpha=.12)
        axes[row, 0].set_ylabel("Group range (km)")
    axes[0, 1].legend(loc="upper right")
    axes[0, 0].set_ylim(700, max(2200, math.ceil(panel_maximum / 100.0) * 100.0))
    for ax in axes[-1]:
        ax.set_xlabel("Frequency (MHz)")
    ionogram_figure = FIGURES / "oblique_wave_600km_ionograms.png"
    fig.savefig(ionogram_figure, dpi=220)
    plt.close(fig)

    altitude, truth, starting, fitted, true_peak = sections["link_midpoint"]
    latitudes = midpoint_points[:, 0]
    true_peaks = np.max(truth, axis=1)
    start_peaks = np.max(starting, axis=1)
    fit_peaks = np.max(fitted, axis=1)
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True,
                             constrained_layout=True)
    axes[0].plot(latitudes, true_peaks / 1e5, "k-", label="IRI-2016 wave truth")
    axes[0].plot(latitudes, start_peaks / 1e5, "o--", color="#777777",
                 label="Starting vertical retrieval")
    axes[0].plot(latitudes, fit_peaks / 1e5, "o-", color="#257e56",
                 label="Oblique-selected retrieval")
    axes[0].set_ylabel("F2 peak density ($10^5$ cm$^{-3}$)")
    axes[0].legend()
    axes[1].plot(latitudes, 100 * (start_peaks / true_peaks - 1), "o--",
                 color="#777777", label="Starting")
    axes[1].plot(latitudes, 100 * (fit_peaks / true_peaks - 1), "o-",
                 color="#257e56", label="Oblique-selected")
    axes[1].axhline(0, color="k", linewidth=.8)
    axes[1].set(xlabel="Link midpoint latitude (degrees)",
                ylabel="Peak-density error (%)")
    for ax in axes:
        ax.grid(alpha=.2)
    peak_figure = FIGURES / "oblique_wave_600km_peaks.png"
    fig.savefig(peak_figure, dpi=220)
    plt.close(fig)

    mask = (altitude >= 150) & (altitude <= 600)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharex=True, sharey=True,
                             constrained_layout=True)
    levels = np.linspace(0, max(np.max(truth[:, mask]), np.max(fitted[:, mask])) / 1e5, 23)
    for ax, values, title in ((axes[0], truth, "Independent density truth"),
                              (axes[1], fitted, "Oblique-selected retrieval")):
        image = ax.contourf(latitudes, altitude, values.T / 1e5,
                            levels=levels, cmap="viridis", extend="max")
        ax.set(title=title, xlabel="Link midpoint latitude (degrees)", ylim=(150, 600))
    fig.colorbar(image, ax=axes[:2], label="Electron density ($10^5$ cm$^{-3}$)",
                 shrink=.8, pad=.01)
    residual = 100 * (fitted - truth) / true_peak[:, None]
    limit = max(10, np.ceil(np.max(abs(residual[:, mask])) / 5) * 5)
    image = axes[2].contourf(latitudes, altitude, residual.T,
                             levels=np.linspace(-limit, limit, 21),
                             cmap="RdBu_r", extend="both")
    axes[2].set(title="Retrieved minus truth", xlabel="Link midpoint latitude (degrees)")
    fig.colorbar(image, ax=axes[2], label="Difference (% of local truth peak)",
                 shrink=.8, pad=.01)
    axes[0].set_ylabel("Altitude (km)")
    density_figure = FIGURES / "oblique_wave_600km_density.png"
    fig.savefig(density_figure, dpi=220)
    plt.close(fig)

    mid = metrics["link_midpoint"]
    tx = metrics["transmitter"]
    if name == "baseline":
        correction_summary = "The starting density was retained."
    else:
        candidate = json.loads((RUN / name / "candidate.json").read_text())
        height_limits = candidate["height_shift_range_km"]
        peak_limits = candidate["correction_range_percent"]
        peak_description = (
            "It applies no additional peak-density scaling."
            if candidate["gain"] == 0.0 else
            f"Its peak-density correction spans {peak_limits[0]:+.1f}% to "
            f"{peak_limits[1]:+.1f}%."
        )
        correction_summary = (
            f"The selected update shifts local profiles by {height_limits[0]:+.1f} "
            f"to {height_limits[1]:+.1f} km. {peak_description}"
        )
    alternate_summary = ""
    if other_candidates:
        alternative = min(other_candidates,
                          key=lambda row: row["combined_observable_score"])
        alternate_summary = (
            f"The almost tied **{alternative['name']}** update scored "
            f"{alternative['combined_observable_score']:.4f} on the observables. "
            f"Its post-selection midpoint peak error was "
            f"{alternative['link_midpoint_peak_density_mae_percent']:.2f}% and "
            f"its profile RMS error was "
            f"{100 * alternative['link_midpoint_profile_mean_nrmse_150_to_600_km']:.2f}%. "
            "These truth-density numbers were not used to choose the retrieval."
        )
    report = f"""# Two-satellite oblique retrieval of the latitude wave

## Setup

Twenty northbound links use the saved 20° latitude wave case. Both spacecraft are at 800 km altitude, with **600 km physical along-track separation** and a nominal northward speed of 8 km/s. The frozen ionosphere is traced in O and X modes from 2 to 10 MHz at 0.1 MHz spacing. The receiver homing gate is 1 km. A path must descend at least 100 km below the spacecraft before it can count as an ionospheric return; direct satellite-to-satellite rays are excluded. All {result['truth_accepted_returns']} accepted O/X returns are saved. Their launch elevations span {result['truth_launch_elevation_range_deg'][0]:.1f}° to {result['truth_launch_elevation_range_deg'][1]:.1f}° and their modeled spacecraft Doppler spans **{result['truth_doppler_range_hz'][0]:.1f} to {result['truth_doppler_range_hz'][1]:.1f} Hz**.

## Retrieval

The starting density is the previous vertical-ionogram retrieval. Smooth latitude corrections were estimated from oblique O/X cutoff differences and 2.5–4.5 MHz group-range residuals. The latter propose a local reflector-height shift. [A six-link screen](data/oblique_wave_600km/pilot_selection.json) compared six updates using the O/X ionogram score alone; its two best candidates were then forward-traced across all 20 links. [The frozen full-pass selection](data/oblique_wave_600km/selection.json) chose **{name}** using 80% O/X ionogram score and 20% Doppler score after 1 Hz quantization. Truth density was not read during either selection.

{correction_summary}

| Mean score across 20 links (lower is better) | Starting | Selected |
| --- | ---: | ---: |
| O/X ionogram | {baseline['mean_ionogram_score']:.4f} | {chosen['mean_ionogram_score']:.4f} |
| 1 Hz Doppler | {baseline['mean_1hz_doppler_score']:.4f} | {chosen['mean_1hz_doppler_score']:.4f} |
| Combined | {baseline['combined_score']:.4f} | {chosen['combined_score']:.4f} |

The ionogram score weights symmetric accepted-return distance (55%), mode-resolved return cutoff (25%), and common-frequency median group-range residual (20%).

The next figure shows four truth/retrieved pairs. O is blue, X red; all accepted paths are shown at their traced 0.1 MHz frequencies and rounded to 1 km group-range bins.

![Truth and retrieved oblique O/X ionograms](figures/oblique_wave_600km_ionograms.png)

## Density check after selection

| Peak-density error across 20 positions | Starting | Selected |
| --- | ---: | ---: |
| Link midpoints, mean absolute | {mid['baseline_peak_density_mae_percent']:.2f}% | {mid['retrieved_peak_density_mae_percent']:.2f}% |
| Transmitter positions, mean absolute | {tx['baseline_peak_density_mae_percent']:.2f}% | {tx['retrieved_peak_density_mae_percent']:.2f}% |
| Link midpoints, 150–600 km normalized profile RMS | {100*mid['baseline_profile_mean_nrmse_150_to_600_km']:.2f}% | {100*mid['retrieved_profile_mean_nrmse_150_to_600_km']:.2f}% |

The selected oblique update improves the overall density profile but increases F2 peak-density error. {alternate_summary}

![Peak densities and errors at the link midpoints](figures/oblique_wave_600km_peaks.png)

![Truth, retrieval, and density residual along the links](figures/oblique_wave_600km_density.png)

## Provenance and limits

The electron-density truth is a separate Fortran IRI-2016 grid with an imposed 900 km latitude wave. The starting retrieval is a transformed PyIRI field. Truth and modeled ionograms share the PyLap ray tracer and homing procedure. The [evaluation record](data/oblique_wave_600km/evaluation.json) was generated after selection. This is an exploratory augmentation of a vertical retrieval that had already been developed on this wave case, so it is not a blind validation. Doppler uses each satellite's local tangent velocity through a static ionosphere. Returns are matched by mode, frequency, and nearest group range within 40 km; unmatched returns incur the maximum Doppler cost, and a 20 Hz difference reaches that maximum. The 1 Hz check is quantization, not a measured noise model.
"""
    (ROOT / "reports/oblique_wave_600km.md").write_text(report)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    single = commands.add_parser("trace")
    single.add_argument("--grid", type=Path, required=True)
    single.add_argument("--index", type=int, required=True)
    single.add_argument("--output", type=Path, required=True)
    single.add_argument("--density-source", required=True)
    many = commands.add_parser("batch")
    many.add_argument("--grid", type=Path, required=True)
    many.add_argument("--output-dir", type=Path, required=True)
    many.add_argument("--workers", type=int, default=1)
    many.add_argument("--indices", nargs="+", type=int, default=list(range(1, 21)))
    many.add_argument("--density-source", required=True)
    candidate = commands.add_parser("build")
    candidate.add_argument("--gain", type=float, required=True)
    candidate.add_argument("--vertical-sigma-km", type=float, required=True)
    candidate.add_argument("--height-gain", type=float, default=0.0)
    candidate.add_argument("--height-fit", choices=("basis", "spline"),
                           default="basis")
    selection = commands.add_parser("select")
    selection.add_argument("--candidates", nargs="+", required=True)
    pilot = commands.add_parser("screen")
    pilot.add_argument("--candidates", nargs="+", required=True)
    pilot.add_argument("--indices", nargs="+", type=int, required=True)
    commands.add_parser("evaluate")
    args = parser.parse_args()
    if args.command == "trace":
        with local_ray_lock():
            result = trace(args.grid, args.index, args.output, args.density_source)
    elif args.command == "batch":
        result = batch(args.grid, args.output_dir, args.workers, args.indices,
                       args.density_source)
    elif args.command == "build":
        result = build_candidate(args.gain, args.vertical_sigma_km,
                                 args.height_gain, args.height_fit)
    elif args.command == "select":
        result = select(args.candidates, list(range(1, 21)))
    elif args.command == "screen":
        result = select(args.candidates, args.indices, RUN / "pilot_selection.json",
                        doppler_weight=0.0)
    else:
        result = evaluate()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
