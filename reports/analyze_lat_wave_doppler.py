"""Assess 8 km/s spacecraft Doppler as an O/X ionogram constraint."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment

from doppler_returns import unique_returns

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reports/data"
FIGURES = ROOT / "reports/figures"
TRUTH = DATA / "lat_wave_doppler_truth_ionograms"
FIT = DATA / "lat_wave_doppler_fit_ionograms"
OUTPUT = DATA / "lat_wave_doppler_assessment.json"


def read(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as source:
        records = np.asarray(source["records"], dtype=float)
        doppler = np.asarray(source["spacecraft_doppler_hz"], dtype=float)
        launch = np.asarray(source["launch_angles_deg"], dtype=float)
        arrival = np.asarray(source["arrival_angles_deg"], dtype=float)
        if (len(records) != len(doppler) or launch.shape != (len(records), 2)
                or arrival.shape != (len(records), 2)):
            raise ValueError(f"Doppler arrays do not align with returns in {path}")
        if (float(source["spacecraft_speed_mps"]) != 8000
                or float(source["spacecraft_track_bearing_deg"]) != 0):
            raise ValueError(f"Unexpected spacecraft kinematics in {path}")
        return records, doppler, launch, arrival


def analyze() -> None:
    all_truth = []
    all_residuals = []
    all_range_residuals = []
    all_frequencies = []
    per_profile = []
    truth_by_profile = {}
    fit_by_profile = {}
    raw_truth_count = 0
    raw_fit_count = 0
    unique_fit_count = 0
    for index in range(1, 21):
        truth = read(TRUTH / f"ionogram_{index:02d}.npz")
        fit = read(FIT / f"ionogram_{index:02d}.npz")
        truth_by_profile[index] = truth
        fit_by_profile[index] = fit
        raw_truth_count += len(truth[0])
        raw_fit_count += len(fit[0])
        truth_unique = unique_returns(truth[0], truth[1])
        fit_unique = unique_returns(fit[0], fit[1])
        all_truth.extend(truth_unique[:, 3])
        unique_fit_count += len(fit_unique)
        matching = []
        for mode in (1, -1):
            for frequency_index in range(81):
                ta = np.flatnonzero((truth_unique[:, 1] == mode)
                                    & (truth_unique[:, 0] == frequency_index))
                fa = np.flatnonzero((fit_unique[:, 1] == mode)
                                    & (fit_unique[:, 0] == frequency_index))
                if not len(ta) or not len(fa):
                    continue
                difference = abs(truth_unique[ta, 2, None]
                                 - fit_unique[None, fa, 2])
                rows, cols = linear_sum_assignment(difference)
                for row, col in zip(rows, cols):
                    if difference[row, col] <= 30.0:
                        ti, fi = ta[row], fa[col]
                        matching.append((float(truth_unique[ti, 3] - fit_unique[fi, 3]),
                                         float(truth_unique[ti, 2] - fit_unique[fi, 2]),
                                         float(2 + 0.1 * frequency_index), mode))
        if matching:
            row = np.asarray(matching)
            all_residuals.extend(row[:, 0])
            all_range_residuals.extend(row[:, 1])
            all_frequencies.extend(row[:, 2])
            per_profile.append({
                "index": index,
                "matched_returns_within_30_km": len(row),
                "median_absolute_doppler_difference_hz": float(np.median(abs(row[:, 0]))),
                "median_absolute_group_range_difference_km": float(np.median(abs(row[:, 1]))),
            })
    signal = np.asarray(all_truth)
    residuals = np.asarray(all_residuals)
    range_residuals = np.asarray(all_range_residuals)
    within_10_km = abs(range_residuals) <= 10.0
    discrepant_1_hz = abs(residuals[within_10_km]) > 1.0
    result = {
        "spacecraft_speed_mps": 8000.0,
        "spacecraft_track_bearing_deg": 0.0,
        "ionosphere_temporally_frozen": True,
        "truth_return_count": raw_truth_count,
        "retrieved_return_count": raw_fit_count,
        "truth_distinct_path_count_after_near_duplicate_grouping": len(signal),
        "retrieved_distinct_path_count_after_near_duplicate_grouping": unique_fit_count,
        "matched_return_count_within_30_km": len(residuals),
        "truth_absolute_doppler_median_hz": float(np.median(abs(signal))),
        "truth_absolute_doppler_p90_hz": float(np.quantile(abs(signal), .9)),
        "truth_absolute_doppler_max_hz": float(np.max(abs(signal))),
        "matched_doppler_difference_median_absolute_hz": float(np.median(abs(residuals))),
        "matched_doppler_difference_p90_absolute_hz": float(np.quantile(abs(residuals), .9)),
        "matched_doppler_difference_rms_hz": float(np.sqrt(np.mean(residuals ** 2))),
        "matched_group_range_difference_median_absolute_km": float(np.median(abs(range_residuals))),
        "range_doppler_residual_correlation": float(np.corrcoef(range_residuals, residuals)[0, 1]),
        "matched_returns_within_10_km": int(np.count_nonzero(within_10_km)),
        "within_10_km_but_doppler_difference_above_1_hz": int(np.count_nonzero(discrepant_1_hz)),
        "within_10_km_but_doppler_difference_above_1_hz_fraction": float(
            np.mean(discrepant_1_hz)),
        "near_duplicate_grouping_range_tolerance_km": 1.0,
        "near_duplicate_grouping_doppler_tolerance_hz": 0.1,
        "per_profile": per_profile,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True,
                             sharex=True, sharey=True)
    for row, index in enumerate((10, 14)):
        for col, (name, data) in enumerate((("Truth", truth_by_profile[index]),
                                             ("Retrieved", fit_by_profile[index]))):
            records, doppler, _, _ = data
            for mode, color, label in ((1, "#135b9a", "O"), (-1, "#b63836", "X")):
                selected = records[:, 1] == mode
                axes[row, col].scatter(2.0 + .1 * records[selected, 0],
                                       doppler[selected], s=11, color=color,
                                       alpha=.75, label=label)
            axes[row, col].set(title=f"{name}, profile {index:02d}", xlim=(2, 10))
            axes[row, col].grid(alpha=.2)
        axes[row, 0].set_ylabel("Spacecraft Doppler (Hz)")
    for axis in axes[-1]:
        axis.set_xlabel("Frequency (MHz)")
    axes[0, 1].legend()
    fig.savefig(FIGURES / "lat_wave_doppler_examples.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    axes[0].hist(signal, bins=np.arange(np.floor(signal.min()) - 1,
                                       np.ceil(signal.max()) + 2, 1), color="#345c81")
    axes[0].set(xlabel="Truth spacecraft Doppler (Hz)", ylabel="Distinct paths",
                title="8 km/s northbound pass, frozen ionosphere")
    axes[1].scatter(range_residuals, residuals, s=8, alpha=.25, color="#69508a")
    axes[1].axhspan(-1, 1, color="#85ad96", alpha=.2, label="±1 Hz gate")
    axes[1].set(xlabel="Truth minus retrieved group range (km)",
                ylabel="Truth minus retrieved Doppler (Hz)",
                title="Matched O/X returns within 30 km")
    axes[1].legend()
    fig.savefig(FIGURES / "lat_wave_doppler_diagnostics.png", dpi=220)
    plt.close(fig)
    print(json.dumps({k: v for k, v in result.items() if k != "per_profile"}, indent=2))


if __name__ == "__main__":
    analyze()
