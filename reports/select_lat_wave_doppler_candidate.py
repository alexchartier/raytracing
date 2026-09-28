"""Rank latitude-wave forward candidates using O/X ionograms and Doppler only.

The candidate set is specified before running this script. No density truth is
loaded here. A 1 Hz Doppler standard deviation is assumed for this experiment.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from doppler_returns import unique_returns
from fit_d_ionogram import Ionogram, score

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reports/data"
TRUTH_IONOGRAMS = DATA / "lat_wave_doppler_truth_ionograms"
BASE_IONOGRAMS = DATA / "lat_wave_doppler_fit_ionograms"
CANDIDATES = DATA / "lat_wave_doppler_candidates"
OUTPUT = DATA / "lat_wave_doppler_selection.json"
ASSUMED_DOPPLER_SIGMA_HZ = 1.0
RANGE_MATCH_GATE_KM = 30.0


def load(path: Path) -> tuple[Ionogram, np.ndarray]:
    ionogram = Ionogram.read(path)
    with np.load(path, allow_pickle=False) as source:
        doppler = np.asarray(source["spacecraft_doppler_hz"], dtype=float)
        speed = float(source["spacecraft_speed_mps"])
        bearing = float(source["spacecraft_track_bearing_deg"])
    if len(doppler) != len(ionogram.records) or not np.all(np.isfinite(doppler)):
        raise ValueError(f"Missing or misaligned Doppler in {path}")
    if speed != 8000.0 or bearing != 0.0:
        raise ValueError(f"Unexpected kinematics in {path}")
    return ionogram, doppler


def doppler_score(observed: tuple[Ionogram, np.ndarray],
                  modeled: tuple[Ionogram, np.ndarray]) -> dict[str, float]:
    """Pair same-frequency O/X returns by range, then penalize Doppler mismatch.

    Every unmatched accepted return contributes the maximum clipped penalty.
    Doppler residuals are capped at 3 assumed standard deviations so one bad
    homing branch cannot dominate a profile.
    """
    obs, obs_doppler = observed
    pred, pred_doppler = modeled
    obs_rows = unique_returns(obs.records, obs_doppler)
    pred_rows = unique_returns(pred.records, pred_doppler)
    penalties = []
    matched_residuals = []
    matched_ranges = []
    matched_count = 0
    for mode in (1, -1):
        for frequency_index in range(81):
            oi = np.flatnonzero((obs_rows[:, 1] == mode)
                                & (obs_rows[:, 0] == frequency_index))
            pi = np.flatnonzero((pred_rows[:, 1] == mode)
                                & (pred_rows[:, 0] == frequency_index))
            if not len(oi) and not len(pi):
                continue
            denominator = max(len(oi), len(pi))
            if not len(oi) or not len(pi):
                penalties.extend([1.0] * denominator)
                continue
            range_difference = abs(obs_rows[oi, 2, None]
                                   - pred_rows[None, pi, 2])
            rows, cols = linear_sum_assignment(range_difference)
            good = range_difference[rows, cols] <= RANGE_MATCH_GATE_KM
            matched_count += int(np.count_nonzero(good))
            differences = abs(obs_rows[oi[rows[good]], 3]
                              - pred_rows[pi[cols[good]], 3])
            matched_residuals.extend(differences.tolist())
            matched_ranges.extend(range_difference[rows[good], cols[good]].tolist())
            penalties.extend(np.minimum(
                differences / (3.0 * ASSUMED_DOPPLER_SIGMA_HZ), 1.0).tolist())
            penalties.extend([1.0] * (denominator - int(np.count_nonzero(good))))
    return {
        "doppler_cost": float(np.mean(penalties)),
        "range_matched_returns": matched_count,
        "unique_observed_returns": len(obs_rows),
        "unique_modeled_returns": len(pred_rows),
        "matched_doppler_mae_hz": (float(np.mean(matched_residuals))
                                   if matched_residuals else None),
        "matched_doppler_median_absolute_hz": (float(np.median(matched_residuals))
                                               if matched_residuals else None),
        "matched_range_mae_km": (float(np.mean(matched_ranges))
                                 if matched_ranges else None),
    }


def main() -> None:
    variant_list = json.loads((CANDIDATES / "candidates.json").read_text())
    variants = [{"name": "baseline", "ionograms": str(BASE_IONOGRAMS.relative_to(ROOT))}]
    variants.extend({"name": variant["name"],
                     "ionograms": str((CANDIDATES / variant["name"] / "recovered").relative_to(ROOT)),
                     "parameter": variant["parameter"], "offset": variant["offset"]}
                    for variant in variant_list)
    observed = [load(TRUTH_IONOGRAMS / f"ionogram_{index:02d}.npz")
                for index in range(1, 21)]
    results = []
    for variant in variants:
        folder = ROOT / variant["ionograms"]
        modeled = [load(folder / f"ionogram_{index:02d}.npz")
                   for index in range(1, 21)]
        per_profile = []
        for index, (truth, fit) in enumerate(zip(observed, modeled), 1):
            ionogram = score(truth[0], fit[0])
            doppler = doppler_score(truth, fit)
            per_profile.append({"index": index, "ionogram_score": ionogram["total"],
                                **doppler})
        row = dict(variant)
        row["mean_ionogram_score"] = float(np.mean([p["ionogram_score"] for p in per_profile]))
        row["mean_doppler_cost"] = float(np.mean([p["doppler_cost"] for p in per_profile]))
        row["matched_return_count"] = int(sum(p["range_matched_returns"] for p in per_profile))
        row["per_profile"] = per_profile
        results.append(row)
    weights = (0.0, 0.1, 0.2, 0.4)
    choices = {}
    for weight in weights:
        key = f"doppler_weight_{weight:g}"
        choices[key] = min(results, key=lambda row: row["mean_ionogram_score"]
                           + weight * row["mean_doppler_cost"])["name"]
    result = {
        "selection_uses_truth_density": False,
        "observations": str(TRUTH_IONOGRAMS.relative_to(ROOT)),
        "assumed_doppler_sigma_hz": ASSUMED_DOPPLER_SIGMA_HZ,
        "range_match_gate_km": RANGE_MATCH_GATE_KM,
        "doppler_cost_cap_sigma": 3.0,
        "primary_doppler_weight": 0.2,
        "selected": choices["doppler_weight_0.2"],
        "choices_by_weight": choices,
        "candidates": results,
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"selected": result["selected"],
                      "choices_by_weight": choices,
                      "scores": [{k: row[k] for k in ("name", "mean_ionogram_score",
                                                  "mean_doppler_cost", "matched_return_count")}
                                 for row in results]}, indent=2))


if __name__ == "__main__":
    main()
