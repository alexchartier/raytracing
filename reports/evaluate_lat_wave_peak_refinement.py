"""Evaluate an ionogram-selected peak correction against withheld wave density.

This is a post-selection check. The selector never opens the truth density.
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

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "reports") not in sys.path:
    sys.path.insert(0, str(ROOT / "reports"))
from evaluate_lat_wave_retrieval import density_at_positions, draw_ionogram  # noqa: E402
from fit_d_ionogram import Ionogram, score  # noqa: E402

DATA = ROOT / "reports/data"
FIGURES = ROOT / "reports/figures"
MANIFEST = DATA / "lat_wave_pass_manifest.json"
TRUTH_DENSITY = DATA / "lat_wave_pass_truth_density.npz"
BASE_DENSITY = DATA / "lat_wave_pass_retrieved_density.npz"
TRUTH_IONOGRAMS = DATA / "lat_wave_final_truth_ionograms"
BASE_IONOGRAMS = DATA / "lat_wave_final_previous_fit_ionograms"
CANDIDATES = DATA / "lat_wave_peak_candidates"


def summarize(selection_path: Path, output_path: Path, figure_prefix: str) -> dict:
    selection = json.loads(selection_path.read_text())
    if selection["selection_uses_truth_density"] or selection["candidate_stage"] != "final":
        raise ValueError("Expected a final ionogram-only candidate selection")
    chosen = selection["chosen_candidate"]
    if chosen == "baseline":
        fit_density, fit_ionograms = BASE_DENSITY, BASE_IONOGRAMS
    else:
        fit_density = CANDIDATES / f"{chosen}_density.npz"
        fit_ionograms = CANDIDATES / chosen / "final"
    manifest = json.loads(MANIFEST.read_text())
    profiles = manifest["profiles"]
    points = np.array([[p["latitude_deg"], p["longitude_deg"]] for p in profiles])
    latitude = points[:, 0]
    altitude, truth = density_at_positions(TRUTH_DENSITY, points)
    base_altitude, base = density_at_positions(BASE_DENSITY, points)
    fit_altitude, retrieved = density_at_positions(fit_density, points)
    if not np.array_equal(altitude, base_altitude) or not np.array_equal(altitude, fit_altitude):
        raise ValueError("Density altitude grids differ")
    truth_peak = np.max(truth, axis=1)
    base_peak = np.max(base, axis=1)
    fit_peak = np.max(retrieved, axis=1)
    base_error = 100 * (base_peak / truth_peak - 1)
    fit_error = 100 * (fit_peak / truth_peak - 1)
    candidate_peak_mae = {}
    for candidate in json.loads((CANDIDATES / "candidates.json").read_text())["candidates"]:
        _, candidate_density = density_at_positions(ROOT / candidate["density"], points)
        candidate_error = 100 * (np.max(candidate_density, axis=1) / truth_peak - 1)
        candidate_peak_mae[candidate["name"]] = float(np.mean(abs(candidate_error)))
    mask = (altitude >= 150) & (altitude <= 600)
    bottomside = (altitude >= 150) & (altitude <= 240)
    topside = (altitude >= 260) & (altitude <= 600)
    base_nrmse = np.sqrt(np.mean((base[:, mask] - truth[:, mask]) ** 2, axis=1)) / truth_peak
    fit_nrmse = np.sqrt(np.mean((retrieved[:, mask] - truth[:, mask]) ** 2, axis=1)) / truth_peak

    # Recompute these metrics from the selected final ionograms to catch a
    # mismatch between the saved selection and the density being evaluated.
    scores = []
    noses = []
    for profile in profiles:
        filename = f"ionogram_{profile['index']:02d}.npz"
        observed = Ionogram.read(TRUTH_IONOGRAMS / filename)
        modeled = Ionogram.read(fit_ionograms / filename)
        scores.append(score(observed, modeled)["total"])
        noses.append([abs(observed.nose(mode) - modeled.nose(mode)) for mode in (1, -1)])
    selected_row = next(row for row in selection["candidates"] if row["name"] == chosen)
    if abs(np.mean(scores) - selected_row["mean_score"]) > 1e-10:
        raise ValueError("Selected ionogram scores differ from saved selection")
    result = {
        "selection_uses_truth_density": False,
        "truth_density_read_only_after_selection": True,
        "selection_file": str(selection_path.resolve().relative_to(ROOT)),
        "chosen_candidate": chosen,
        "truth_density": str(TRUTH_DENSITY.relative_to(ROOT)),
        "retrieved_density": str(fit_density.relative_to(ROOT)),
        "truth_ionograms": str(TRUTH_IONOGRAMS.relative_to(ROOT)),
        "retrieved_ionograms": str(fit_ionograms.relative_to(ROOT)),
        "profile_count": len(profiles),
        "baseline_peak_density_mae_percent": float(np.mean(abs(base_error))),
        "selected_peak_density_mae_percent": float(np.mean(abs(fit_error))),
        "selected_peak_density_max_absolute_error_percent": float(np.max(abs(fit_error))),
        "selected_peak_density_median_absolute_error_percent": float(np.median(abs(fit_error))),
        "profiles_with_smaller_peak_error_than_baseline": int(np.count_nonzero(abs(fit_error) < abs(base_error))),
        "post_selection_peak_density_mae_percent_by_candidate": candidate_peak_mae,
        "baseline_peak_density_error_percent_by_profile": base_error.tolist(),
        "selected_peak_density_error_percent_by_profile": fit_error.tolist(),
        "baseline_fof2_mae_mhz": float(np.mean(abs(0.00898 * (np.sqrt(base_peak) - np.sqrt(truth_peak))))),
        "selected_fof2_mae_mhz": float(np.mean(abs(0.00898 * (np.sqrt(fit_peak) - np.sqrt(truth_peak))))),
        "baseline_density_profile_mean_nrmse_150_to_600_km": float(np.mean(base_nrmse)),
        "selected_density_profile_mean_nrmse_150_to_600_km": float(np.mean(fit_nrmse)),
        "selected_bottomside_mean_nrmse_150_to_240_km": float(np.mean(
            np.sqrt(np.mean((retrieved[:, bottomside] - truth[:, bottomside]) ** 2,
                            axis=1)) / truth_peak)),
        "selected_topside_mean_nrmse_260_to_600_km": float(np.mean(
            np.sqrt(np.mean((retrieved[:, topside] - truth[:, topside]) ** 2,
                            axis=1)) / truth_peak)),
        "selected_ionogram_mean_score": float(np.mean(scores)),
        "selected_nose_mae_mhz_by_mode": np.mean(noses, axis=0).tolist(),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2) + "\n")
    FIGURES.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True, sharex=True)
    axes[0].plot(latitude, truth_peak / 1e5, "k-", label="Independent truth")
    axes[0].plot(latitude, base_peak / 1e5, "o--", color="#777777", label="Previous retrieval")
    axes[0].plot(latitude, fit_peak / 1e5, "o-", color="#237a55", label="Selected peak correction")
    axes[0].set(ylabel="Peak density (100,000 cm$^{-3}$)", title="F2 peak across the 20-position wave pass")
    axes[0].legend()
    axes[1].plot(latitude, base_error, "o--", color="#777777", label="Previous retrieval")
    axes[1].plot(latitude, fit_error, "o-", color="#237a55", label="Selected peak correction")
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set(xlabel="Latitude (degrees)", ylabel="Peak density error (%)")
    axes[1].legend()
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.savefig(FIGURES / f"{figure_prefix}_peaks.png", dpi=220)
    plt.close(fig)

    indices = (3, 10, 14, 18)
    fig, axes = plt.subplots(len(indices), 2, figsize=(13, 13), constrained_layout=True,
                             sharex=True, sharey=True)
    for row, index in enumerate(indices):
        filename = f"ionogram_{index:02d}.npz"
        draw_ionogram(axes[row, 0], TRUTH_IONOGRAMS / filename,
                      f"Truth, profile {index:02d} ({latitude[index - 1]:.1f}°)")
        draw_ionogram(axes[row, 1], fit_ionograms / filename,
                      f"Retrieved, profile {index:02d}")
        axes[row, 0].set_ylabel("Group range (km)")
    for ax in axes[-1]:
        ax.set_xlabel("Frequency (MHz)")
    axes[0, 1].legend(loc="upper right", frameon=True)
    fig.savefig(FIGURES / f"{figure_prefix}_ionograms.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True,
                             sharex=True, sharey=True)
    levels = np.linspace(0, max(np.max(truth[:, mask]), np.max(retrieved[:, mask])) / 1e5, 23)
    for ax, values, title in ((axes[0], truth, "Independent truth"),
                              (axes[1], retrieved, "Selected retrieval")):
        image = ax.contourf(latitude, altitude, values.T / 1e5, levels=levels,
                            cmap="viridis", extend="max")
        ax.set(title=title, xlabel="Latitude (degrees)", ylim=(150, 600))
    fig.colorbar(image, ax=axes[:2], label="Electron density (100,000 cm$^{-3}$)",
                 shrink=0.8, pad=0.01)
    difference = 100 * (retrieved - truth) / truth_peak[:, None]
    limit = max(10, np.ceil(np.max(abs(difference[:, mask])) / 5) * 5)
    residual = axes[2].contourf(latitude, altitude, difference.T,
                                 levels=np.linspace(-limit, limit, 21),
                                 cmap="RdBu_r", extend="both")
    fig.colorbar(residual, ax=axes[2], label="Difference (% of local truth peak)",
                 shrink=0.8, pad=0.01)
    axes[2].set(title="Retrieved minus truth", xlabel="Latitude (degrees)")
    axes[0].set_ylabel("Altitude (km)")
    fig.savefig(FIGURES / f"{figure_prefix}_density.png", dpi=220)
    plt.close(fig)
    print(json.dumps({key: value for key, value in result.items()
                      if not key.endswith("_by_profile")}, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("selection", type=Path)
    parser.add_argument("--output", type=Path, default=CANDIDATES / "evaluation.json")
    parser.add_argument("--figure-prefix", default="lat_wave_peak_refinement")
    arguments = parser.parse_args()
    summarize(arguments.selection, arguments.output, arguments.figure_prefix)
