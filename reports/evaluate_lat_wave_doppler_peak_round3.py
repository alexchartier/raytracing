"""Check a frozen low-Doppler ionogram selection against truth density."""

from __future__ import annotations

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
from select_lat_wave_doppler_peak_round3 import gated_nose  # noqa: E402

DATA = ROOT / "reports/data"
ROUND = DATA / "lat_wave_doppler_peak_round3"
FIGURES = ROOT / "reports/figures"
SELECTION = ROUND / "final_selection.json"
TRUTH_DENSITY = DATA / "lat_wave_pass_truth_density.npz"
BASE_DENSITY = DATA / "lat_wave_joint_round2/peak_full_density.npz"
TRUTH_IONOGRAMS = DATA / "lat_wave_final_truth_ionograms"
BASE_IONOGRAMS = DATA / "lat_wave_joint_round2/peak_full/final"


def quantized_gated_nose(path: Path, mode: int) -> float:
    with np.load(path, allow_pickle=False) as data:
        records = np.asarray(data["records"], dtype=float)
        frequencies = np.asarray(data["frequencies_mhz"], dtype=float)
        doppler = np.rint(np.asarray(data["spacecraft_doppler_hz"], dtype=float))
    selected = records[(records[:, 1] == mode) & (abs(doppler) <= 15.0)]
    if not len(selected):
        raise ValueError(f"No 1 Hz quantized gated return in {path}")
    return float(frequencies[int(np.max(selected[:, 0]))])


def evaluate() -> dict:
    selection = json.loads(SELECTION.read_text())
    if (selection["selection_uses_truth_density"]
            or selection["stage"] != "final"
            or selection["profile_indices"] != list(range(1, 21))):
        raise ValueError("Expected frozen ionogram-only selection over all 20 profiles")
    name = selection["chosen_candidate"]
    density_path = (BASE_DENSITY if name == "baseline" else
                    ROUND / f"{name}_density.npz")
    ionograms = (BASE_IONOGRAMS if name == "baseline" else ROUND / name / "final")
    profiles = json.loads((DATA / "lat_wave_pass_manifest.json").read_text())["profiles"]
    points = np.array([[row["latitude_deg"], row["longitude_deg"]]
                       for row in profiles])
    latitude = points[:, 0]
    altitude, truth = density_at_positions(TRUTH_DENSITY, points)
    base_altitude, baseline = density_at_positions(BASE_DENSITY, points)
    fit_altitude, fit = density_at_positions(density_path, points)
    if not (np.array_equal(altitude, base_altitude)
            and np.array_equal(altitude, fit_altitude)):
        raise ValueError("Density altitude grids differ")
    true_peak, base_peak, fit_peak = (np.max(array, axis=1)
                                      for array in (truth, baseline, fit))
    base_error = 100 * (base_peak / true_peak - 1)
    fit_error = 100 * (fit_peak / true_peak - 1)
    mask = (altitude >= 150) & (altitude <= 600)
    base_nrmse = np.sqrt(np.mean((baseline[:, mask] - truth[:, mask]) ** 2, axis=1)) / true_peak
    fit_nrmse = np.sqrt(np.mean((fit[:, mask] - truth[:, mask]) ** 2, axis=1)) / true_peak
    selected_row = next(row for row in selection["candidates"] if row["name"] == name)
    baseline_row = next(row for row in selection["candidates"] if row["name"] == "baseline")
    score_change = np.array([a["total"] - b["total"] for a, b in zip(
        selected_row["scores_by_profile"], baseline_row["scores_by_profile"])])
    rng = np.random.default_rng(20260928)
    bootstrap = rng.choice(score_change, size=(20_000, 20), replace=True).mean(axis=1)
    quantized_gate_error = {}
    for candidate, directory in (("baseline", BASE_IONOGRAMS), ("selected", ionograms)):
        quantized_gate_error[candidate] = float(np.mean([
            abs(quantized_gated_nose(TRUTH_IONOGRAMS / f"ionogram_{i:02d}.npz", mode)
                - quantized_gated_nose(directory / f"ionogram_{i:02d}.npz", mode))
            for i in range(1, 21) for mode in (1, -1)]))
    result = {
        "selection_uses_truth_density": False,
        "method_development_used_prior_truth_evaluation": True,
        "current_candidate_selection_used_ionograms_only": True,
        "selected_candidate": name,
        "candidate_selection": str(SELECTION.relative_to(ROOT)),
        "retrieved_density": str(density_path.relative_to(ROOT)),
        "baseline_peak_density_mae_percent": float(np.mean(abs(base_error))),
        "selected_peak_density_mae_percent": float(np.mean(abs(fit_error))),
        "baseline_peak_density_max_absolute_error_percent": float(np.max(abs(base_error))),
        "selected_peak_density_max_absolute_error_percent": float(np.max(abs(fit_error))),
        "selected_peak_density_median_absolute_error_percent": float(np.median(abs(fit_error))),
        "baseline_fof2_mae_mhz": float(np.mean(abs(.00898 *
                                                   (np.sqrt(base_peak) - np.sqrt(true_peak))))),
        "selected_fof2_mae_mhz": float(np.mean(abs(.00898 *
                                                   (np.sqrt(fit_peak) - np.sqrt(true_peak))))),
        "baseline_profile_mean_nrmse_150_to_600_km": float(np.mean(base_nrmse)),
        "selected_profile_mean_nrmse_150_to_600_km": float(np.mean(fit_nrmse)),
        "baseline_standard_ionogram_score": baseline_row["standard_score"],
        "selected_standard_ionogram_score": selected_row["standard_score"],
        "baseline_combined_ionogram_score": baseline_row["combined_score"],
        "selected_combined_ionogram_score": selected_row["combined_score"],
        "baseline_gated_cutoff_mae_mhz_by_mode": baseline_row["gated_cutoff_mae_mhz_by_mode"],
        "selected_gated_cutoff_mae_mhz_by_mode": selected_row["gated_cutoff_mae_mhz_by_mode"],
        "baseline_gated_cutoff_mae_with_1hz_doppler_quantization_mhz":
            quantized_gate_error["baseline"],
        "selected_gated_cutoff_mae_with_1hz_doppler_quantization_mhz":
            quantized_gate_error["selected"],
        "standard_score_change_bootstrap_95_percent_interval": np.quantile(
            bootstrap, [.025, .975]).tolist(),
        "profiles_with_smaller_peak_error": int(np.count_nonzero(abs(fit_error) < abs(base_error))),
        "baseline_peak_error_percent_by_profile": base_error.tolist(),
        "selected_peak_error_percent_by_profile": fit_error.tolist(),
    }
    (ROUND / "evaluation.json").write_text(json.dumps(result, indent=2) + "\n")
    FIGURES.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
    axes[0].plot(latitude, true_peak / 1e5, "k-", label="Independent synthetic truth")
    axes[0].plot(latitude, base_peak / 1e5, "o--", color="#777777", label="Previous retrieval")
    axes[0].plot(latitude, fit_peak / 1e5, "o-", color="#257e56", label="Selected update")
    axes[0].set(ylabel="Peak density ($10^5$ cm$^{-3}$)",
                title="F2 peak across the 20-position pass")
    axes[0].legend()
    axes[1].plot(latitude, base_error, "o--", color="#777777", label="Previous retrieval")
    axes[1].plot(latitude, fit_error, "o-", color="#257e56", label="Selected update")
    axes[1].axhline(0, color="black", linewidth=.8)
    axes[1].set(xlabel="Latitude (degrees)", ylabel="Peak-density error (%)")
    axes[1].legend()
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(FIGURES / "lat_wave_doppler_peak_round3_peaks.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(4, 2, figsize=(13, 13), sharex=True, sharey=True,
                             constrained_layout=True)
    for row, index in enumerate((3, 8, 14, 18)):
        filename = f"ionogram_{index:02d}.npz"
        draw_ionogram(axes[row, 0], TRUTH_IONOGRAMS / filename,
                      f"Truth, profile {index:02d} ({latitude[index - 1]:.1f}°)")
        draw_ionogram(axes[row, 1], ionograms / filename,
                      f"Retrieved, profile {index:02d}")
        axes[row, 0].set_ylabel("Group range (km)")
    for ax in axes[-1]:
        ax.set_xlabel("Frequency (MHz)")
    axes[0, 1].legend(loc="upper right")
    fig.savefig(FIGURES / "lat_wave_doppler_peak_round3_ionograms.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharex=True, sharey=True,
                             constrained_layout=True)
    levels = np.linspace(0, max(np.max(truth[:, mask]), np.max(fit[:, mask])) / 1e5, 23)
    for ax, values, title in ((axes[0], truth, "Independent synthetic truth"),
                              (axes[1], fit, "Selected retrieval")):
        image = ax.contourf(latitude, altitude, values.T / 1e5, levels=levels,
                            cmap="viridis", extend="max")
        ax.set(title=title, xlabel="Latitude (degrees)", ylim=(150, 600))
    fig.colorbar(image, ax=axes[:2], label="Electron density ($10^5$ cm$^{-3}$)",
                 shrink=.8, pad=.01)
    residual = 100 * (fit - truth) / true_peak[:, None]
    limit = max(10, np.ceil(np.max(abs(residual[:, mask])) / 5) * 5)
    image = axes[2].contourf(latitude, altitude, residual.T,
                             levels=np.linspace(-limit, limit, 21),
                             cmap="RdBu_r", extend="both")
    fig.colorbar(image, ax=axes[2], label="Difference (% of local truth peak)",
                 shrink=.8, pad=.01)
    axes[2].set(title="Retrieved minus truth", xlabel="Latitude (degrees)")
    axes[0].set_ylabel("Altitude (km)")
    fig.savefig(FIGURES / "lat_wave_doppler_peak_round3_density.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True, constrained_layout=True)
    for row, mode in enumerate((1, -1)):
        truth_nose = [gated_nose(TRUTH_IONOGRAMS / f"ionogram_{i:02d}.npz", mode)
                      for i in range(1, 21)]
        base_nose = [gated_nose(BASE_IONOGRAMS / f"ionogram_{i:02d}.npz", mode)
                     for i in range(1, 21)]
        fit_nose = [gated_nose(ionograms / f"ionogram_{i:02d}.npz", mode)
                    for i in range(1, 21)]
        axes[row].plot(latitude, truth_nose, "k-", label="Observed")
        axes[row].plot(latitude, base_nose, "o--", color="#777777", label="Previous")
        axes[row].plot(latitude, fit_nose, "o-", color="#257e56", label="Selected")
        axes[row].set(ylabel=f"{('O', 'X')[row]} cutoff (MHz)")
        axes[row].grid(alpha=.2)
    axes[0].set_title("Highest accepted return with |spacecraft Doppler| ≤ 15 Hz")
    axes[0].legend()
    axes[-1].set_xlabel("Latitude (degrees)")
    fig.savefig(FIGURES / "lat_wave_doppler_peak_round3_gated_cutoffs.png", dpi=220)
    plt.close(fig)
    print(json.dumps({key: value for key, value in result.items()
                      if not key.endswith("_by_profile")}, indent=2))
    return result


if __name__ == "__main__":
    evaluate()
