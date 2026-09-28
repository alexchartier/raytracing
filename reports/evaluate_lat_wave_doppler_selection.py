"""Evaluate the already selected Doppler candidate against withheld density."""

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

DATA = ROOT / "reports/data"
FIGURES = ROOT / "reports/figures"
SELECTION = DATA / "lat_wave_doppler_selection.json"
MANIFEST = DATA / "lat_wave_pass_manifest.json"
TRUTH = DATA / "lat_wave_pass_truth_density.npz"
BASELINE = DATA / "lat_wave_pass_retrieved_density.npz"
CANDIDATES = DATA / "lat_wave_doppler_candidates"
OUTPUT = DATA / "lat_wave_doppler_retrieval_evaluation.json"


def metrics(density: np.ndarray, truth: np.ndarray) -> dict[str, float]:
    peak = np.max(density, axis=1)
    truth_peak = np.max(truth, axis=1)
    error = 100 * (peak / truth_peak - 1)
    return {
        "peak_density_mean_absolute_error_percent": float(np.mean(abs(error))),
        "peak_density_median_absolute_error_percent": float(np.median(abs(error))),
        "fof2_mean_absolute_error_mhz": float(np.mean(abs(
            0.00898 * np.sqrt(peak) - 0.00898 * np.sqrt(truth_peak)))),
        "peak_density_signed_error_percent_by_profile": error.tolist(),
    }


def main() -> None:
    selection = json.loads(SELECTION.read_text())
    chosen = selection["selected"]
    manifest = json.loads(MANIFEST.read_text())
    points = np.array([[p["latitude_deg"], p["longitude_deg"]]
                       for p in manifest["profiles"]])
    latitudes = points[:, 0]
    selected_path = (BASELINE if chosen == "baseline"
                     else CANDIDATES / f"{chosen}_density.npz")
    altitude, truth = density_at_positions(TRUTH, points)
    base_altitude, baseline = density_at_positions(BASELINE, points)
    selected_altitude, selected = density_at_positions(selected_path, points)
    if not (np.array_equal(altitude, base_altitude)
            and np.array_equal(altitude, selected_altitude)):
        raise ValueError("Altitude axes differ")
    result = {
        "selection_fixed_before_truth_density_loaded": True,
        "selection_file": str(SELECTION.relative_to(ROOT)),
        "selected": chosen,
        "baseline": metrics(baseline, truth),
        "selected_candidate": metrics(selected, truth),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    truth_peak = np.max(truth, axis=1)
    baseline_peak = np.max(baseline, axis=1)
    selected_peak = np.max(selected, axis=1)
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                             constrained_layout=True)
    axes[0].plot(latitudes, truth_peak / 1e5, "k-", label="Withheld IRI truth")
    axes[0].plot(latitudes, baseline_peak / 1e5, color="#135b9a",
                 label="Ionogram-only baseline")
    axes[0].plot(latitudes, selected_peak / 1e5, color="#b63836",
                 label=f"Best nearby candidate: {chosen}")
    axes[0].set(ylabel="Peak density (100,000 cm$^{-3}$)",
                title="Latitude-wave pass: nearby-candidate comparison")
    axes[0].legend()
    axes[1].plot(latitudes, 100 * (baseline_peak / truth_peak - 1),
                 "o-", color="#135b9a", label="Ionogram-only baseline")
    axes[1].plot(latitudes, 100 * (selected_peak / truth_peak - 1),
                 "o-", color="#b63836", label="Best nearby candidate")
    axes[1].axhline(0, color="k", linewidth=.7)
    axes[1].set(xlabel="Latitude (°)", ylabel="Peak density error (%)")
    for axis in axes:
        axis.grid(alpha=.2)
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES / "lat_wave_doppler_retrieval_peaks.png", dpi=220)
    plt.close(fig)

    truth_ionograms = DATA / "lat_wave_doppler_truth_ionograms"
    selected_ionograms = (DATA / "lat_wave_doppler_fit_ionograms" if chosen == "baseline"
                          else CANDIDATES / chosen / "recovered")
    indices = (3, 10, 14, 18)
    fig, axes = plt.subplots(len(indices), 2, figsize=(13, 13), constrained_layout=True,
                             sharex=True, sharey=True)
    for row, index in enumerate(indices):
        draw_ionogram(axes[row, 0],
                      truth_ionograms / f"ionogram_{index:02d}.npz",
                      f"Truth, profile {index:02d} ({latitudes[index - 1]:.1f}°)")
        draw_ionogram(axes[row, 1],
                      selected_ionograms / f"ionogram_{index:02d}.npz",
                      f"Selected, profile {index:02d}")
        axes[row, 0].set_ylabel("Group range (km)")
    for axis in axes[-1]:
        axis.set_xlabel("Frequency (MHz)")
    axes[0, 1].legend(loc="upper right", frameon=True)
    fig.savefig(FIGURES / "lat_wave_doppler_retrieval_ionograms.png", dpi=220)
    plt.close(fig)
    print(json.dumps({k: v for k, v in result.items() if k not in
                      ("baseline", "selected_candidate")}, indent=2))
    print(json.dumps({"baseline": {k: v for k, v in result["baseline"].items()
                                   if not k.endswith("by_profile")},
                      "selected_candidate": {
                          k: v for k, v in result["selected_candidate"].items()
                          if not k.endswith("by_profile")}}, indent=2))


if __name__ == "__main__":
    main()
