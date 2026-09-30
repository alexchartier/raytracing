"""Post-selection density and ionogram check for general-field IRI pilots.

The full-ray score and candidate choice must be frozen before this command
opens the independent IRI-2016 wave density.
"""

from __future__ import annotations

import argparse
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
from python_raytrace.grid import load_ionosphere_grid_netcdf  # noqa: E402

DATA = ROOT / "reports/data"
RUN = DATA / "general_field_iri_wave_structured"
TRUTH = DATA / "lat_wave_pass_truth_density.npz"
PREVIOUS = DATA / "lat_wave_doppler_peak_round3/wide_density.npz"
OBSERVED = DATA / "lat_wave_final_truth_ionograms"
PREVIOUS_IONOGRAMS = DATA / "lat_wave_doppler_peak_round3/wide/final"


def _profiles(path: Path, locations: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if path.suffix == ".nc":
        grid = load_ionosphere_grid_netcdf(path)
        latitude, longitude, altitude, density = (
            grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km,
            grid.iono_en_grid)
    else:
        with np.load(path, allow_pickle=False) as source:
            latitude = source["latitudes_deg"]
            longitude = source["longitudes_deg"]
            altitude = source["altitudes_km"]
            density = source["electron_density_cm3"]
    sampled = RegularGridInterpolator((latitude, longitude), density,
                                      bounds_error=True)(locations)
    return np.asarray(altitude), np.asarray(sampled)


def _draw_ionogram(ax, path: Path, title: str) -> None:
    with np.load(path, allow_pickle=False) as source:
        records = np.asarray(source["records"], dtype=float)
        frequencies = np.asarray(source["frequencies_mhz"], dtype=float)
    for mode, color, label in ((1, "#135b9a", "O"), (-1, "#b63836", "X")):
        selected = records[records[:, 1] == mode]
        ax.scatter(frequencies[selected[:, 0].astype(int)],
                   np.round(selected[:, 2]), s=7, marker="s", linewidths=0,
                   color=color, alpha=.78, label=label)
    ax.set(xlim=(2, 10), ylim=(150, 2100), title=title)
    ax.grid(alpha=.12)


def evaluate(selection_path: Path, score_path: Path, name: str | None) -> None:
    selection = json.loads(selection_path.read_text())
    scores = json.loads(score_path.read_text())
    if (selection["selection_uses_truth_density"]
            or scores["selection_uses_truth_density"]
            or sorted(selection["profile_indices"]) != list(range(1, 21))
            or selection["profile_indices"] != scores["profile_indices"]
            or Path(selection["full_ray_scores"]).resolve() != score_path.resolve()):
        raise ValueError("Expected a frozen all-profile ionogram-only selection")
    scored = {row["name"]: row for row in scores["scores"]}
    chosen = min(scores["scores"], key=lambda row: row["mean_combined_score"])
    if selection["chosen_candidate"] != chosen["name"]:
        raise ValueError("Frozen choice does not match the ionogram scores")
    if name is None:
        name = selection["chosen_candidate"]
    plan = json.loads((RUN / "plan.json").read_text())
    if name == "previous":
        grid_path = PREVIOUS
        modeled_ionograms = PREVIOUS_IONOGRAMS
    elif name == "combined":
        grid_path = RUN / "combined/grid.nc"
        modeled_ionograms = RUN / "forward/combined/ionograms"
    else:
        candidate = next((row for row in plan["candidates"] if row["name"] == name), None)
        if candidate is None:
            raise ValueError(f"Unknown candidate {name}")
        grid_path = ROOT / candidate["grid"]
        modeled_ionograms = RUN / "forward" / name / "ionograms"
    score_row = scored[name]
    previous_score = scored["previous"]
    manifest = json.loads((DATA / "lat_wave_pass_manifest.json").read_text())["profiles"]
    locations = np.array([(row["latitude_deg"], row["longitude_deg"])
                          for row in manifest])
    latitude = locations[:, 0]
    altitude, truth = _profiles(TRUTH, locations)
    previous_altitude, previous = _profiles(PREVIOUS, locations)
    baseline_altitude, baseline = _profiles(RUN / "baseline/grid.nc", locations)
    fit_altitude, retrieved = _profiles(grid_path, locations)
    if not (np.array_equal(altitude, previous_altitude)
            and np.array_equal(altitude, baseline_altitude)
            and np.array_equal(altitude, fit_altitude)):
        raise ValueError("Density altitude grids differ")
    true_peak = np.max(truth, axis=1)
    true_height = altitude[np.argmax(truth, axis=1)]
    mask = (altitude >= 150) & (altitude <= 600)
    middlemask = (altitude >= 220) & (altitude <= 600)
    topmask = (altitude >= 220) & (altitude <= 800)

    def metrics(profiles: np.ndarray) -> dict:
        peak = profiles.max(axis=1)
        height = altitude[np.argmax(profiles, axis=1)]
        local = np.array([np.interp(800.0, altitude, row) for row in profiles])
        local_truth = np.array([np.interp(800.0, altitude, row) for row in truth])
        peak_error = 100 * (peak / true_peak - 1)
        height_error = height - true_height
        lower_error = (100 * np.sqrt(np.mean((profiles[:, mask] - truth[:, mask]) ** 2,
                                             axis=1)) / true_peak)
        middle_error = (100 * np.sqrt(np.mean(
            (profiles[:, middlemask] - truth[:, middlemask]) ** 2,
            axis=1)) / true_peak)
        upper_error = (100 * np.sqrt(np.mean((profiles[:, topmask] - truth[:, topmask]) ** 2,
                                             axis=1)) / true_peak)
        return {
            "peak_density_mae_percent": float(np.mean(abs(peak_error))),
            "peak_density_bias_percent": float(np.mean(peak_error)),
            "fof2_mae_mhz": float(np.mean(abs(.00898 * (np.sqrt(peak) - np.sqrt(true_peak))))),
            "hmf2_mae_km": float(np.mean(abs(height_error))),
            "hmf2_bias_km": float(np.mean(height_error)),
            "density_nrmse_150_to_600_percent": float(np.mean(lower_error)),
            "density_nrmse_220_to_600_percent": float(np.mean(middle_error)),
            "density_nrmse_220_to_800_percent": float(np.mean(upper_error)),
            "local_800km_density_mae_percent": float(np.mean(abs(
                100 * (local / local_truth - 1)))),
            "peak_density_error_percent_by_profile": peak_error.tolist(),
            "hmf2_error_km_by_profile": height_error.tolist(),
            "density_nrmse_150_to_600_percent_by_profile": lower_error.tolist(),
            "density_nrmse_220_to_600_percent_by_profile": middle_error.tolist(),
            "density_nrmse_220_to_800_percent_by_profile": upper_error.tolist(),
        }

    result = {
        "candidate": name,
        "selection_uses_truth_density": False,
        "prior_method_development_used_truth_diagnostics": True,
        "observed_density_source": "independent Fortran IRI-2016 with imposed latitude wave",
        "local_800km_density_is_an_assimilated_observation": True,
        "candidate_grid": str(grid_path.relative_to(ROOT)),
        "candidate_ionograms": str(modeled_ionograms.relative_to(ROOT)),
        "previous": metrics(previous),
        "candidate_metrics": metrics(retrieved),
        "previous_combined_ionogram_score": previous_score["mean_combined_score"],
        "candidate_combined_ionogram_score": score_row["mean_combined_score"],
    }
    output = RUN / f"evaluation_{name}.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    output.chmod(0o600)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for row, index in enumerate((8, 14)):
        ax = axes[row, 0]
        ax.plot(truth[index - 1] / 1e5, altitude, "k", label="IRI wave truth")
        ax.plot(previous[index - 1] / 1e5, altitude, "--", color="#777777",
                label="Previous retrieval")
        ax.plot(retrieved[index - 1] / 1e5, altitude, color="#257e56",
                label="New candidate")
        ax.set(xlim=(0, 5), ylim=(150, 850),
               xlabel="Density ($10^5$ cm$^{-3}$)", ylabel="Altitude (km)",
               title=f"Profile {index:02d}")
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8)
    for ax, key, title, ylabel in (
            (axes[0, 1], "peak_density_error_percent_by_profile", "Peak density", "Error (%)"),
            (axes[1, 1], "hmf2_error_km_by_profile", "Peak height", "Error (km)")):
        ax.plot(latitude, result["previous"][key], "o--", color="#777777")
        ax.plot(latitude, result["candidate_metrics"][key], "o-", color="#257e56")
        ax.axhline(0, color="black", linewidth=.8)
        ax.set(title=title, xlabel="Latitude (degrees)", ylabel=ylabel)
        ax.grid(alpha=.2)
    fig.savefig(RUN / f"density_comparison_{name}.png", dpi=190)
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.5), sharex=True,
                             sharey=True, constrained_layout=True)
    density_levels = np.linspace(0.0, 4.6, 24)
    for ax, field, title in (
            (axes[0, 0], truth, "Independent IRI wave truth"),
            (axes[0, 1], retrieved, f"Retrieved: {name}")):
        density_image = ax.contourf(latitude, altitude, (field / 1e5).T,
                                    levels=density_levels, cmap="viridis",
                                    extend="max")
        ax.set_title(title)
    fig.colorbar(density_image, ax=axes[0, :],
                 label="Electron density ($10^5$ cm$^{-3}$)")
    residual = (retrieved - truth) / 1e5
    residual_image = axes[1, 0].contourf(
        latitude, altitude, residual.T, levels=np.linspace(-0.8, 0.8, 25),
        cmap="RdBu_r", extend="both")
    axes[1, 0].set_title("Retrieved − truth")
    fig.colorbar(residual_image, ax=axes[1, 0],
                 label="Density difference ($10^5$ cm$^{-3}$)")
    increment = 100 * (retrieved / baseline - 1)
    increment_image = axes[1, 1].contourf(
        latitude, altitude, increment.T, levels=np.linspace(0, 3, 21),
        cmap="magma", extend="both")
    axes[1, 1].set_title(f"{name} − anchored baseline")
    fig.colorbar(increment_image, ax=axes[1, 1],
                 label="Density change (%)")
    for ax in axes.flat:
        ax.set(xlim=(latitude.min(), latitude.max()), ylim=(150, 800))
    for ax in axes[:, 0]:
        ax.set_ylabel("Altitude (km)")
    for ax in axes[1, :]:
        ax.set_xlabel("Latitude (degrees)")
    fig.savefig(RUN / f"latitude_altitude_density_{name}.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(4, 2, figsize=(12, 12), sharex=True, sharey=True,
                             constrained_layout=True)
    for row, index in enumerate((3, 8, 14, 18)):
        filename = f"ionogram_{index:02d}.npz"
        _draw_ionogram(axes[row, 0], OBSERVED / filename,
                       f"Truth {index:02d} ({latitude[index - 1]:.1f}°)")
        _draw_ionogram(axes[row, 1], modeled_ionograms / filename,
                       f"Retrieved {index:02d}")
        axes[row, 0].set_ylabel("Group range (km)")
    for ax in axes[-1]:
        ax.set_xlabel("Frequency (MHz)")
    axes[0, 1].legend(loc="upper right")
    fig.savefig(RUN / f"ionogram_comparison_{name}.png", dpi=190)
    plt.close(fig)
    print(json.dumps({key: value for key, value in result.items()
                      if key not in ("previous", "candidate_metrics")}, indent=2))
    print(json.dumps({"previous": {k: v for k, v in result["previous"].items()
                                    if not k.endswith("_by_profile")},
                      "candidate": {k: v for k, v in result["candidate_metrics"].items()
                                    if not k.endswith("_by_profile")}}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--candidate-name")
    args = parser.parse_args()
    evaluate(args.selection, args.scores, args.candidate_name)
