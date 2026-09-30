"""Evaluate a frozen ionogram-only SAMI retrieval against withheld density truth."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "reports"))
from python_raytrace.grid import load_ionosphere_grid_netcdf  # noqa: E402
from ionogram_metrics import Ionogram, score  # noqa: E402

CASE = ROOT / "reports/data/sami3_wave_20170111_0600"
RUN = ROOT / "reports/data/general_field_sami3_wave"


def profiles(path: Path, locations: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    grid = load_ionosphere_grid_netcdf(path)
    return np.asarray(grid.altitudes_km), np.asarray(RegularGridInterpolator(
        (grid.latitudes_deg, grid.longitudes_deg), grid.iono_en_grid,
        bounds_error=True)(locations))


def metrics(altitude: np.ndarray, truth: np.ndarray, model: np.ndarray) -> dict:
    peak_truth = truth.max(axis=1)
    peak_model = model.max(axis=1)
    true_height = altitude[np.argmax(truth, axis=1)]
    model_height = altitude[np.argmax(model, axis=1)]
    mask = (altitude >= 150) & (altitude <= 800)
    upper = (altitude >= 350) & (altitude <= 800)
    peak_error = 100 * (peak_model / peak_truth - 1)
    h_error = model_height - true_height
    rmse = np.sqrt(np.mean((model[:, mask] - truth[:, mask]) ** 2, axis=1))
    upper_rmse = np.sqrt(np.mean((model[:, upper] - truth[:, upper]) ** 2, axis=1))
    return {
        "peak_density_mae_percent": float(np.mean(abs(peak_error))),
        "peak_density_bias_percent": float(np.mean(peak_error)),
        "fof2_mae_mhz": float(np.mean(abs(.00898 * (np.sqrt(peak_model) - np.sqrt(peak_truth))))),
        "hmf2_mae_km": float(np.mean(abs(h_error))),
        "hmf2_bias_km": float(np.mean(h_error)),
        "profile_nrmse_150_to_800_percent": float(np.mean(100 * rmse / peak_truth)),
        "topside_nrmse_350_to_800_percent": float(np.mean(100 * upper_rmse / peak_truth)),
        "peak_density_error_percent_by_profile": peak_error.tolist(),
        "hmf2_error_km_by_profile": h_error.tolist(),
    }


def wave_metrics(latitudes: np.ndarray, truth: np.ndarray, model: np.ndarray) -> dict:
    """Compare along-track foF2 structure after removing a linear background."""
    frequencies = [.00898 * np.sqrt(field.max(axis=1)) for field in (truth, model)]
    design = np.column_stack((np.ones(len(latitudes)), latitudes - latitudes.mean()))
    residuals = [value - design @ np.linalg.lstsq(design, value, rcond=None)[0]
                 for value in frequencies]
    actual, fitted = residuals
    true_rms = float(np.sqrt(np.mean(actual ** 2)))
    model_rms = float(np.sqrt(np.mean(fitted ** 2)))
    return {
        "fof2_residual_correlation": float(np.corrcoef(actual, fitted)[0, 1]),
        "truth_residual_rms_mhz": true_rms,
        "retrieved_residual_rms_mhz": model_rms,
        "retrieved_to_truth_rms_ratio": model_rms / true_rms,
        "residual_mae_mhz": float(np.mean(abs(fitted - actual))),
    }


def draw_ionogram(ax, path: Path, title: str) -> None:
    with np.load(path, allow_pickle=False) as item:
        frequency = item["frequencies_mhz"]
        records = item["records"]
    for mode, color, label in ((1, "#2166ac", "O"), (-1, "#b2182b", "X")):
        selected = records[records[:, 1] == mode]
        ax.scatter(frequency[selected[:, 0].astype(int)], selected[:, 2],
                   s=10, c=color, marker="s", linewidths=0, label=label)
    ax.set(xlim=(2, 10), ylim=(950, 2050), xlabel="Frequency (MHz)",
           ylabel="Group range (km)", title=title)
    ax.grid(alpha=.2)


def main() -> None:
    selection = json.loads((RUN / "selection.json").read_text())
    scores = json.loads((RUN / "full_scores.json").read_text())
    plan = json.loads((RUN / "plan.json").read_text())
    indices = [1, 5, 9, 12, 14, 17, 20]
    if (selection["selection_uses_truth_density"] or
            selection["profile_indices"] != indices or
            scores["profile_indices"] != indices or
            selection["chosen_candidate"] != min(
                scores["scores"], key=lambda row: row["mean_combined_score"])["name"]):
        raise ValueError("Expected the frozen seven-ionogram selection")
    name = selection["chosen_candidate"]
    candidate = next((row for row in plan["candidates"] if row["name"] == name), None)
    grid_path = (CASE / "pyiri_prior_grid.nc" if name == "prior" else
                 RUN / "combined/grid.nc" if name == "combined" else
                 ROOT / candidate["grid"])
    original_rows = json.loads((CASE / "manifest.json").read_text())["profiles"]
    locations = np.array([(row["latitude_deg"], row["longitude_deg"])
                          for row in original_rows])
    altitude, truth = profiles(CASE / "truth_grid.nc", locations)
    prior_altitude, prior = profiles(CASE / "pyiri_prior_grid.nc", locations)
    retrieved_altitude, retrieved = profiles(grid_path, locations)
    if not (np.array_equal(altitude, prior_altitude)
            and np.array_equal(altitude, retrieved_altitude)):
        raise ValueError("Density altitude axes differ")
    observed_rows = np.array([i - 1 for i in indices])
    result = {
        "selection_uses_truth_density": False,
        "truth_model": "SAMI3/HIAMCM remapped snapshot",
        "prior_model": "independent PyIRI",
        "local_density_observations": "seven synthetic in-situ 800 km samples",
        "observed_profile_indices": indices,
        "selected_candidate": name,
        "selected_grid": str(grid_path.relative_to(ROOT)),
        "selected_ionogram_score": selection["chosen_mean_combined_score"],
        "prior_ionogram_score": float(np.mean([
            score(Ionogram.read(CASE / "vertical_truth" / f"ionogram_{i:02d}.npz"),
                  Ionogram.read(CASE / "vertical_prior" / f"ionogram_{i:02d}.npz"))["total"]
            for i in indices])),
        "observed_positions": {
            "prior": metrics(altitude, truth[observed_rows], prior[observed_rows]),
            "retrieved": metrics(altitude, truth[observed_rows], retrieved[observed_rows]),
        },
        "full_twenty_position_track": {
            "prior": metrics(altitude, truth, prior),
            "retrieved": metrics(altitude, truth, retrieved),
            "prior_wave_structure": wave_metrics(locations[:, 0], truth, prior),
            "retrieved_wave_structure": wave_metrics(locations[:, 0], truth, retrieved),
        },
    }
    output = RUN / "evaluation.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    output.chmod(0o600)

    latitudes = locations[:, 0]
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), constrained_layout=True,
                             sharex=True, sharey=True)
    scale = (4.0, 5.7)
    for ax, field, title in zip(axes[:2], (truth, retrieved),
                                ("SAMI3/HIAMCM truth", "Retrieved from PyIRI")):
        mesh = ax.contourf(latitudes, altitude, np.log10(np.maximum(field.T, 1.0)),
                           levels=np.linspace(*scale, 18), cmap="viridis", extend="both")
        ax.scatter(latitudes[observed_rows], np.full(len(indices), 790),
                   s=12, color="white", edgecolors="black", linewidths=.3)
        ax.set(title=title, xlabel="Latitude (degrees)", ylim=(150, 800))
    fig.colorbar(mesh, ax=axes[:2], label="log10 electron density (cm$^{-3}$)")
    difference = 100 * (retrieved - truth) / truth.max(axis=1)[:, None]
    diff = axes[2].contourf(latitudes, altitude, difference.T,
                            levels=np.linspace(-40, 40, 17), cmap="RdBu_r", extend="both")
    axes[2].set(title="Difference / true peak (%)", xlabel="Latitude (degrees)")
    fig.colorbar(diff, ax=axes[2], label="Percent")
    axes[0].set_ylabel("Altitude (km)")
    fig.savefig(RUN / "density_cut.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    index = 14
    draw_ionogram(axes[0], CASE / "vertical_truth" / f"ionogram_{index:02d}.npz",
                   "SAMI3/HIAMCM truth, profile 14")
    modeled = (CASE / "vertical_prior" if name == "prior" else
               RUN / "forward" / name / "ionograms")
    draw_ionogram(axes[1], modeled / f"ionogram_{index:02d}.npz",
                   "Selected fit, profile 14")
    axes[0].legend()
    fig.savefig(RUN / "ionogram_14.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 5), constrained_layout=True,
                             sharex=True, sharey=True)
    for ax, index in zip(axes, (9, 14)):
        row = index - 1
        ax.plot(truth[row] / 1e5, altitude, color="black", linewidth=2,
                label="SAMI truth")
        ax.plot(prior[row] / 1e5, altitude, color="#777777", linestyle="--",
                label="PyIRI prior")
        ax.plot(retrieved[row] / 1e5, altitude, color="#2b8cbe",
                label="Selected fit")
        ax.set(title=f"Profile {index:02d}", xlabel="Electron density ($10^5$ cm$^{-3}$)",
               xlim=(0, 5), ylim=(150, 800))
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Altitude (km)")
    axes[0].legend()
    fig.savefig(RUN / "density_profiles.png", dpi=180)
    plt.close(fig)
    print(json.dumps({"selected": name,
                      "observed": result["observed_positions"]["retrieved"],
                      "full_track": result["full_twenty_position_track"]["retrieved"]}, indent=2))


if __name__ == "__main__":
    main()
