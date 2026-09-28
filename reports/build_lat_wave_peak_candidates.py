"""Propose smooth 20 km F2 peak corrections from O/X ionogram noses only.

The existing wave density already contains the earlier fixed 3.0646% peak
correction. This script estimates an additional latitude-dependent correction
from the observed-minus-modeled O/X return limits. The fit and selection do
not read the withheld truth density.
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
MANIFEST = DATA / "lat_wave_pass_manifest.json"
FIT = DATA / "lat_wave_pass_initial_fit.json"
OBSERVED = DATA / "lat_wave_final_truth_ionograms"
MODELED = DATA / "lat_wave_final_previous_fit_ionograms"
BASE_GRID = DATA / "lat_wave_pass_retrieved_forward_grid.nc"
OUTPUT = DATA / "lat_wave_peak_candidates"
PEAK_SIGMA_KM = 20.0


def nose_and_support(path: Path, mode: int) -> tuple[float, int]:
    with np.load(path, allow_pickle=False) as source:
        records = np.asarray(source["records"], dtype=float)
        frequencies = np.asarray(source["frequencies_mhz"], dtype=float)
    selected = records[records[:, 1] == mode]
    if not len(selected):
        raise ValueError(f"No mode {mode} returns in {path}")
    last = int(np.max(selected[:, 0]))
    support = int(len(np.unique(selected[selected[:, 0] >= last - 4, 0])))
    return float(frequencies[last]), support


def fit_correction() -> tuple[np.ndarray, dict]:
    manifest = json.loads(MANIFEST.read_text())
    fit = json.loads(FIT.read_text())
    latitudes = np.array([row["latitude_deg"] for row in manifest["profiles"]])
    center = float(np.mean(latitudes))
    halfspan_km = float(np.max(abs(6371.0088 * np.deg2rad(latitudes - center))))
    wavelength_km = float(fit["wave_wavelength_km"])

    def basis(values: np.ndarray) -> np.ndarray:
        x_km = 6371.0088 * np.deg2rad(values - center)
        phase = 2 * np.pi * x_km / wavelength_km
        return np.column_stack((np.ones(len(values)), x_km / halfspan_km,
                                np.cos(phase), np.sin(phase)))

    observed = np.empty((20, 2))
    modeled = np.empty((20, 2))
    support = np.empty((20, 2))
    for index in range(1, 21):
        for column, mode in enumerate((1, -1)):
            observed[index - 1, column], obs_support = nose_and_support(
                OBSERVED / f"ionogram_{index:02d}.npz", mode)
            modeled[index - 1, column], fit_support = nose_and_support(
                MODELED / f"ionogram_{index:02d}.npz", mode)
            support[index - 1, column] = min(1.0, obs_support / 3,
                                              fit_support / 3)
    # For a local peak density multiplier (1+a), plasma frequency changes
    # approximately by (a/2)*f. Forward raytraces will test this first guess.
    design = np.repeat(basis(latitudes), 2, axis=0)
    frequency_sensitivity = 0.5 * modeled.reshape(-1, 1) * design
    residual_mhz = (observed - modeled).reshape(-1)
    regularization = np.array([0.07, 0.07, 0.07, 0.07])

    def objective(coefficients: np.ndarray) -> np.ndarray:
        return np.r_[np.sqrt(support.reshape(-1))
                     * (frequency_sensitivity @ coefficients - residual_mhz),
                     0.05 * coefficients / regularization]

    result = least_squares(objective, np.zeros(4), loss="soft_l1",
                           f_scale=0.12)
    coefficients = np.asarray(result.x, dtype=float)
    diagnostics = {
        "selection_uses_truth_density": False,
        "observed_ionograms": str(OBSERVED.relative_to(ROOT)),
        "baseline_ionograms": str(MODELED.relative_to(ROOT)),
        "baseline_forward_grid": str(BASE_GRID.relative_to(ROOT)),
        "method": "Soft-L1 fit of O/X nose residuals to local plasma-frequency peak sensitivity",
        "correction_basis": "constant + linear latitude + cosine/sine at fitted wave wavelength",
        "latitude_center_deg": center,
        "latitude_halfspan_km": halfspan_km,
        "wave_wavelength_km": wavelength_km,
        "coefficients": coefficients.tolist(),
        "observed_noses_mhz": observed.tolist(),
        "modeled_noses_mhz": modeled.tolist(),
        "nose_support_weights": support.tolist(),
        "estimated_extra_peak_fraction_by_profile":
            (basis(latitudes) @ coefficients).tolist(),
        "peak_sigma_km": PEAK_SIGMA_KM,
    }
    return coefficients, diagnostics


def build() -> None:
    coefficients, diagnostics = fit_correction()
    grid = load_ionosphere_grid_netcdf(BASE_GRID)
    density = np.asarray(grid.iono_en_grid, dtype=float)
    latitudes = np.asarray(grid.latitudes_deg, dtype=float)
    longitudes = np.asarray(grid.longitudes_deg, dtype=float)
    altitudes = np.asarray(grid.altitudes_km, dtype=float)
    center = diagnostics["latitude_center_deg"]
    halfspan_km = diagnostics["latitude_halfspan_km"]
    wavelength_km = diagnostics["wave_wavelength_km"]
    x_km = 6371.0088 * np.deg2rad(latitudes - center)
    phase = 2 * np.pi * x_km / wavelength_km
    latitude_basis = np.column_stack((np.ones(len(latitudes)), x_km / halfspan_km,
                                      np.cos(phase), np.sin(phase)))
    smooth_alpha = latitude_basis @ coefficients
    peak_altitudes = altitudes[np.argmax(density, axis=2)]
    envelope = np.exp(-0.5 * ((altitudes[None, None, :]
                              - peak_altitudes[:, :, None]) / PEAK_SIGMA_KM) ** 2)
    variants = {
        "uniform_3pct": np.full(len(latitudes), 0.03),
        "smooth_half": 0.5 * smooth_alpha,
        "smooth_full": smooth_alpha,
        "smooth_strong": 1.5 * smooth_alpha,
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    candidates = []
    for name, correction in variants.items():
        correction = np.clip(correction, -0.10, 0.10)
        candidate_density = density * (1.0 + correction[:, None, None] * envelope)
        if not np.all(np.isfinite(candidate_density)) or np.any(candidate_density <= 0):
            raise ValueError(f"Invalid candidate density for {name}")
        collision = grid.collision_freq
        if (grid.electron_temp_k is not None and grid.ion_temp_k is not None
                and grid.neutral_species_cm3 is not None):
            collision = effective_collision_frequency(
                grid.electron_temp_k, grid.ion_temp_k,
                candidate_density * 1e6, grid.neutral_species_cm3)
        candidate_grid = replace(grid, iono_en_grid=candidate_density,
                                 iono_en_grid_5=candidate_density,
                                 collision_freq=collision)
        grid_path = OUTPUT / f"{name}_grid.nc"
        density_path = OUTPUT / f"{name}_density.npz"
        save_ionosphere_grid_netcdf(grid_path, candidate_grid)
        np.savez_compressed(density_path, latitudes_deg=latitudes,
                            longitudes_deg=longitudes, altitudes_km=altitudes,
                            electron_density_cm3=candidate_density,
                            model=np.array(f"PyIRI wave fit plus {name} F2 peak correction"))
        candidates.append({"name": name, "grid": str(grid_path.relative_to(ROOT)),
                           "density": str(density_path.relative_to(ROOT)),
                           "peak_sigma_km": PEAK_SIGMA_KM,
                           "min_extra_peak_fraction": float(np.min(correction)),
                           "max_extra_peak_fraction": float(np.max(correction))})
    diagnostics["candidates"] = candidates
    (OUTPUT / "candidates.json").write_text(json.dumps(diagnostics, indent=2) + "\n")
    print(json.dumps({"coefficients": diagnostics["coefficients"],
                      "estimated_extra_peak_fraction_by_profile":
                          diagnostics["estimated_extra_peak_fraction_by_profile"],
                      "candidates": candidates}, indent=2))


if __name__ == "__main__":
    build()
