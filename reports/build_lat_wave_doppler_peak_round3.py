"""Fit a shared peak correction from low-Doppler O/X return limits.

The 15 Hz gate selects near-nadir returns at 8 km/s spacecraft speed. Its
cutoff residuals are smoothed across all 20 profiles in the existing wave
basis, then applied with a conservative quarter step. Candidate generation
does not read electron-density truth.
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from python_raytrace.absorption import effective_collision_frequency  # noqa: E402
from python_raytrace.grid import load_ionosphere_grid_netcdf, save_ionosphere_grid_netcdf  # noqa: E402

DATA = ROOT / "reports/data"
OBSERVED = DATA / "lat_wave_final_truth_ionograms"
MODELED = DATA / "lat_wave_joint_round2/peak_full/final"
START_GRID = DATA / "lat_wave_joint_round2/peak_full_grid.nc"
MANIFEST = DATA / "lat_wave_pass_manifest.json"
WAVE_FIT = DATA / "lat_wave_pass_initial_fit.json"
OUTPUT = DATA / "lat_wave_doppler_peak_round3"
DOPPLER_LIMIT_HZ = 15.0
STEP = .25


def gated_nose(path: Path, mode: int) -> tuple[float, int]:
    with np.load(path, allow_pickle=False) as data:
        records = np.asarray(data["records"], dtype=float)
        frequencies = np.asarray(data["frequencies_mhz"], dtype=float)
        doppler = np.asarray(data["spacecraft_doppler_hz"], dtype=float)
    selected = records[(records[:, 1] == mode) &
                       (abs(doppler) <= DOPPLER_LIMIT_HZ)]
    if not len(selected):
        raise ValueError(f"No low-Doppler mode {mode} return in {path}")
    last = int(np.max(selected[:, 0]))
    support = int(np.count_nonzero(np.unique(selected[:, 0]) >= last - 4))
    return float(frequencies[last]), support


def basis(latitudes: np.ndarray, center: float, halfspan_km: float,
          wavelength_km: float) -> np.ndarray:
    x_km = 6371.0088 * np.deg2rad(latitudes - center)
    phase = 2 * np.pi * x_km / wavelength_km
    return np.column_stack((np.ones(len(latitudes)), x_km / halfspan_km,
                            np.cos(phase), np.sin(phase)))


def fit() -> dict:
    positions = json.loads(MANIFEST.read_text())["profiles"]
    latitudes = np.array([row["latitude_deg"] for row in positions])
    center = float(np.mean(latitudes))
    halfspan_km = float(np.max(abs(6371.0088 * np.deg2rad(latitudes - center))))
    wavelength_km = float(json.loads(WAVE_FIT.read_text())["wave_wavelength_km"])
    design = basis(latitudes, center, halfspan_km, wavelength_km)
    observed = np.empty((20, 2))
    modeled = np.empty((20, 2))
    support = np.empty((20, 2))
    for index in range(1, 21):
        for column, mode in enumerate((1, -1)):
            observed[index - 1, column], a = gated_nose(
                OBSERVED / f"ionogram_{index:02d}.npz", mode)
            modeled[index - 1, column], b = gated_nose(
                MODELED / f"ionogram_{index:02d}.npz", mode)
            support[index - 1, column] = min(a, b, 5) / 5
    proxy = 2 * (observed - modeled) / modeled
    # Average both measured modes, then fit one common spatial correction.
    # Huber-like loss prevents a single dropped near-nadir return dominating.
    target = np.mean(proxy, axis=1)
    weight = np.sqrt(np.mean(support, axis=1))

    def residual(coefficients: np.ndarray) -> np.ndarray:
        return np.r_[weight * (design @ coefficients - target) / .07,
                     coefficients / np.array([.20, .20, .20, .20])]

    solution = least_squares(residual, np.zeros(4), loss="soft_l1", f_scale=1.0)
    correction = design @ solution.x
    return {"selection_uses_truth_density": False,
            "observed_ionograms": str(OBSERVED.relative_to(ROOT)),
            "modeled_ionograms": str(MODELED.relative_to(ROOT)),
            "starting_grid": str(START_GRID.relative_to(ROOT)),
            "doppler_limit_hz": DOPPLER_LIMIT_HZ,
            "spacecraft_speed_mps": 8000.0,
            "latitude_center_deg": center,
            "latitude_halfspan_km": halfspan_km,
            "wave_wavelength_km": wavelength_km,
            "basis": "constant + linear latitude + cosine/sine at fitted wave wavelength",
            "coefficients": solution.x.tolist(),
            "gated_observed_noses_mhz": observed.tolist(),
            "gated_modeled_noses_mhz": modeled.tolist(),
            "gated_nose_support_weights": support.tolist(),
            "unscaled_peak_fraction_by_profile": correction.tolist(),
            "candidate_step": STEP}


def build() -> None:
    result = fit()
    grid = load_ionosphere_grid_netcdf(START_GRID)
    density = np.asarray(grid.iono_en_grid, dtype=float)
    latitudes = np.asarray(grid.latitudes_deg, dtype=float)
    altitudes = np.asarray(grid.altitudes_km, dtype=float)
    design = basis(latitudes, result["latitude_center_deg"],
                   result["latitude_halfspan_km"], result["wave_wavelength_km"])
    fraction = STEP * (design @ np.asarray(result["coefficients"]))
    if np.max(abs(fraction)) > .05:
        raise ValueError("Peak update exceeds 5% trust region")
    peak_altitudes = altitudes[np.argmax(density, axis=2)]
    candidates = []
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for sigma_km in (20.0, 65.0):
        name = "narrow" if sigma_km == 20.0 else "wide"
        envelope = np.exp(-.5 * ((altitudes[None, None, :]
                                 - peak_altitudes[:, :, None]) / sigma_km) ** 2)
        revised = density * (1 + fraction[:, None, None] * envelope)
        if np.any(revised <= 0) or not np.all(np.isfinite(revised)):
            raise ValueError(f"Invalid density for {name}")
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
        np.savez_compressed(density_path, latitudes_deg=grid.latitudes_deg,
                            longitudes_deg=grid.longitudes_deg,
                            altitudes_km=altitudes, electron_density_cm3=revised,
                            model=np.array(f"Low-Doppler joint wave correction, {name}"))
        candidates.append({"name": name, "grid": str(grid_path.relative_to(ROOT)),
                           "density": str(density_path.relative_to(ROOT)),
                           "peak_sigma_km": sigma_km,
                           "min_peak_fraction": float(np.min(fraction)),
                           "max_peak_fraction": float(np.max(fraction))})
    result["candidates"] = candidates
    (OUTPUT / "candidates.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"coefficients": result["coefficients"],
                      "unscaled_peak_fraction_by_profile":
                          result["unscaled_peak_fraction_by_profile"],
                      "candidates": candidates}, indent=2))


if __name__ == "__main__":
    build()
