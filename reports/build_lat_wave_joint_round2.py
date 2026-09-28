"""Build a second joint wave fit from 20 O/X ionogram residuals.

The starting density is the first-round ionogram-selected wavy grid. A shared
latitude wave in peak height is fitted to centered group-range residuals, and
a small peak-density trend is fitted to O/X nose residuals. No truth density is
opened here. Full 3-D ray traces select among the proposed corrections.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "reports") not in sys.path:
    sys.path.insert(0, str(ROOT / "reports"))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from fit_d_ionogram import Ionogram  # noqa: E402
from python_raytrace.absorption import effective_collision_frequency  # noqa: E402
from python_raytrace.grid import load_ionosphere_grid_netcdf, save_ionosphere_grid_netcdf  # noqa: E402

DATA = ROOT / "reports/data"
OBSERVED = DATA / "lat_wave_final_truth_ionograms"
START_IONOGRAMS = DATA / "lat_wave_peak_candidates/smooth_half/final"
START_GRID = DATA / "lat_wave_peak_candidates/smooth_half_grid.nc"
START_DENSITY = DATA / "lat_wave_peak_candidates/smooth_half_density.npz"
MANIFEST = DATA / "lat_wave_pass_manifest.json"
WAVE_FIT = DATA / "lat_wave_pass_initial_fit.json"
OUTPUT = DATA / "lat_wave_joint_round2"
PEAK_SIGMA_KM = 20.0


def basis(latitudes: np.ndarray, wavelength_km: float) -> tuple[np.ndarray, np.ndarray]:
    x_km = 6371.0088 * np.deg2rad(latitudes - np.mean(latitudes))
    x_norm = x_km / max(abs(x_km))
    phase = 2 * np.pi * x_km / wavelength_km
    return (np.column_stack((x_norm, np.cos(phase), np.sin(phase))),
            np.column_stack((np.ones(len(latitudes)), x_norm,
                             np.cos(phase), np.sin(phase))))


def fit_joint_update() -> dict:
    positions = json.loads(MANIFEST.read_text())["profiles"]
    latitudes = np.array([item["latitude_deg"] for item in positions], dtype=float)
    wavelength = float(json.loads(WAVE_FIT.read_text())["wave_wavelength_km"])
    height_basis, peak_basis = basis(latitudes, wavelength)
    observed = [Ionogram.read(OBSERVED / f"ionogram_{index:02d}.npz")
                for index in range(1, 21)]
    modeled = [Ionogram.read(START_IONOGRAMS / f"ionogram_{index:02d}.npz")
               for index in range(1, 21)]

    # Use the latitude variation at a fixed mode/frequency. Removing each
    # mode/frequency median absorbs a systematic forward-model range offset.
    samples: list[tuple[int, int, int, float]] = []
    for position, (truth, fit) in enumerate(zip(observed, modeled)):
        for mode in (1, -1):
            truth_ridge, fit_ridge = truth.ridge(mode), fit.ridge(mode)
            last = min(int(np.max(truth.records[truth.records[:, 1] == mode, 0])),
                       int(np.max(fit.records[fit.records[:, 1] == mode, 0])))
            for index in truth_ridge.keys() & fit_ridge.keys():
                frequency = float(truth.frequencies[index])
                upper = 4.4 if mode == 1 else 4.8
                if 2.5 <= frequency <= upper and index <= last - 4:
                    samples.append((position, mode, index,
                                    truth_ridge[index] - fit_ridge[index]))
    groups: dict[tuple[int, int], list[int]] = defaultdict(list)
    for row, (_, mode, index, _) in enumerate(samples):
        groups[mode, index].append(row)
    retained = sorted(row for members in groups.values() if len(members) >= 8
                      for row in members)
    samples = [samples[row] for row in retained]
    groups = defaultdict(list)
    for row, (_, mode, index, _) in enumerate(samples):
        groups[mode, index].append(row)
    profile_index = np.array([row[0] for row in samples], dtype=int)
    range_residual = np.array([row[3] for row in samples], dtype=float)
    centered_residual = range_residual.copy()
    centered_design = height_basis[profile_index].copy()
    for members in groups.values():
        centered_residual[members] -= np.median(range_residual[members])
        centered_design[members] -= np.median(centered_design[members], axis=0)
    counts = np.bincount(profile_index, minlength=20)
    profile_weight = np.sqrt(np.mean(counts) / np.maximum(counts, 1))

    def range_objective(coefficients: np.ndarray) -> np.ndarray:
        # For a topside reflection, lifting the layer by one km shortens the
        # two-way group path by roughly two km. Ray tracing tests this proxy.
        return np.r_[profile_weight[profile_index]
                     * (-2 * centered_design @ coefficients - centered_residual) / 35.0,
                     coefficients / np.array([25.0, 30.0, 30.0])]

    range_fit = least_squares(range_objective, np.zeros(3),
                              loss="soft_l1", f_scale=1.0)
    height_km = height_basis @ range_fit.x

    nose_residual = np.array([[a.nose(mode) - b.nose(mode) for mode in (1, -1)]
                              for a, b in zip(observed, modeled)])
    modeled_nose = np.array([[item.nose(mode) for mode in (1, -1)]
                             for item in modeled])
    mode_indicators = np.tile(np.eye(2), (20, 1))
    nose_design = np.column_stack((0.5 * modeled_nose.reshape(-1, 1)
                                   * np.repeat(peak_basis, 2, axis=0),
                                   mode_indicators))

    def nose_objective(coefficients: np.ndarray) -> np.ndarray:
        return np.r_[(nose_design @ coefficients - nose_residual.reshape(-1)) / 0.1,
                     coefficients[:4] / np.array([0.03, 0.04, 0.04, 0.04]),
                     coefficients[4:] / 0.15]

    nose_fit = least_squares(nose_objective, np.zeros(6),
                             loss="soft_l1", f_scale=1.0)
    peak_fraction = peak_basis @ nose_fit.x[:4]
    return {
        "selection_uses_truth_density": False,
        "observed_ionograms": str(OBSERVED.relative_to(ROOT)),
        "starting_ionograms": str(START_IONOGRAMS.relative_to(ROOT)),
        "starting_grid": str(START_GRID.relative_to(ROOT)),
        "starting_density": str(START_DENSITY.relative_to(ROOT)),
        "starting_wave_wavelength_km": wavelength,
        "range_samples": len(samples),
        "range_mode_frequency_groups": len(groups),
        "height_basis": "linear latitude + cosine and sine at starting wave wavelength",
        "height_coefficients_km": range_fit.x.tolist(),
        "height_update_km_by_profile": height_km.tolist(),
        "height_proxy": "two-way group-range change approximately -2 km per km layer shift",
        "peak_basis": "constant + linear latitude + cosine and sine at starting wave wavelength",
        "peak_coefficients_fraction": nose_fit.x[:4].tolist(),
        "remaining_mode_nose_offsets_mhz": nose_fit.x[4:].tolist(),
        "peak_update_fraction_by_profile": peak_fraction.tolist(),
        "peak_sigma_km": PEAK_SIGMA_KM,
    }


def build() -> None:
    fit = fit_joint_update()
    grid = load_ionosphere_grid_netcdf(START_GRID)
    density = np.asarray(grid.iono_en_grid, dtype=float)
    latitudes = np.asarray(grid.latitudes_deg, dtype=float)
    longitudes = np.asarray(grid.longitudes_deg, dtype=float)
    altitudes = np.asarray(grid.altitudes_km, dtype=float)
    center = float(np.mean([row["latitude_deg"]
                            for row in json.loads(MANIFEST.read_text())["profiles"]]))
    x_km = 6371.0088 * np.deg2rad(latitudes - center)
    x_norm = x_km / max(abs(6371.0088 * np.deg2rad(
        np.array([row["latitude_deg"]
                  for row in json.loads(MANIFEST.read_text())["profiles"]]) - center)))
    phase = 2 * np.pi * x_km / fit["starting_wave_wavelength_km"]
    height_design = np.column_stack((x_norm, np.cos(phase), np.sin(phase)))
    peak_design = np.column_stack((np.ones(len(latitudes)), x_norm,
                                   np.cos(phase), np.sin(phase)))
    delta_height = height_design @ np.asarray(fit["height_coefficients_km"])
    delta_peak = peak_design @ np.asarray(fit["peak_coefficients_fraction"])
    if np.max(abs(delta_height)) > 20 or np.max(abs(delta_peak)) > 0.05:
        raise ValueError("Joint update exceeds conservative grid bounds")
    candidates = {
        "height_half": (0.5, 0.0),
        "height_full": (1.0, 0.0),
        "height_double": (2.0, 0.0),
        "peak_half": (0.0, 0.5),
        "peak_full": (0.0, 1.0),
        "peak_double": (0.0, 2.0),
        "joint_half": (0.5, 1.0),
        "joint_full": (1.0, 1.0),
        "joint_double": (2.0, 1.0),
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, (height_strength, peak_strength) in candidates.items():
        shifted = np.empty_like(density)
        for i, height in enumerate(delta_height * height_strength):
            for j in range(len(longitudes)):
                source = density[i, j]
                shifted[i, j] = np.interp(altitudes - height, altitudes, source,
                                          left=source[0], right=source[-1])
        local_peak = altitudes[np.argmax(shifted, axis=2)]
        envelope = np.exp(-0.5 * ((altitudes[None, None, :]
                                  - local_peak[:, :, None]) / PEAK_SIGMA_KM) ** 2)
        revised = shifted * (1 + peak_strength * delta_peak[:, None, None] * envelope)
        if np.any(revised <= 0) or not np.all(np.isfinite(revised)):
            raise ValueError(f"Invalid candidate density: {name}")
        collision = grid.collision_freq
        if (grid.electron_temp_k is not None and grid.ion_temp_k is not None
                and grid.neutral_species_cm3 is not None):
            collision = effective_collision_frequency(
                grid.electron_temp_k, grid.ion_temp_k,
                revised * 1e6, grid.neutral_species_cm3)
        candidate_grid = replace(grid, iono_en_grid=revised,
                                 iono_en_grid_5=revised, collision_freq=collision)
        grid_path = OUTPUT / f"{name}_grid.nc"
        density_path = OUTPUT / f"{name}_density.npz"
        save_ionosphere_grid_netcdf(grid_path, candidate_grid)
        np.savez_compressed(density_path, latitudes_deg=latitudes,
                            longitudes_deg=longitudes, altitudes_km=altitudes,
                            electron_density_cm3=revised,
                            model=np.array(f"Second-round joint wave fit: {name}"))
        manifest.append({"name": name, "height_strength": height_strength,
                         "peak_strength": peak_strength,
                         "grid": str(grid_path.relative_to(ROOT)),
                         "density": str(density_path.relative_to(ROOT))})
    fit["candidates"] = manifest
    (OUTPUT / "candidates.json").write_text(json.dumps(fit, indent=2) + "\n")
    print(json.dumps({"height_coefficients_km": fit["height_coefficients_km"],
                      "height_update_km_by_profile": fit["height_update_km_by_profile"],
                      "peak_update_fraction_by_profile": fit["peak_update_fraction_by_profile"],
                      "candidates": manifest}, indent=2))


if __name__ == "__main__":
    build()
