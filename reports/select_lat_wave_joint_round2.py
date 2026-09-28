"""Select a shared-wave second-round candidate from O/X ionograms only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from fit_d_ionogram import Ionogram, score

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reports/data"
ROUND = DATA / "lat_wave_joint_round2"
TRUTH = DATA / "lat_wave_final_truth_ionograms"
RAW_TRUTH = DATA / "lat_wave_pass_ionograms"
BASELINE = DATA / "lat_wave_peak_candidates/smooth_half"


def select(stage: str, names: list[str], indices: list[int], output: Path) -> dict:
    known = {row["name"]: row for row in
             json.loads((ROUND / "candidates.json").read_text())["candidates"]}
    if not names or any(name not in known for name in names):
        raise ValueError("Choose one or more names from candidates.json")
    if not indices or any(index < 1 or index > 20 for index in indices):
        raise ValueError("Profile indices must be between 1 and 20")
    if stage == "final" and sorted(set(indices)) != list(range(1, 21)):
        raise ValueError("Final selection must use the full 20-profile pass")
    truth_dir = RAW_TRUTH if stage == "raw" else TRUTH
    baseline_dir = BASELINE / stage
    directories = {"baseline": baseline_dir}
    directories.update({name: ROUND / name / stage for name in names})
    observed = {index: Ionogram.read(truth_dir / f"ionogram_{index:02d}.npz")
                for index in indices}
    candidates = []
    for name, directory in directories.items():
        scores = []
        nose_errors = []
        for index in indices:
            modeled = Ionogram.read(directory / f"ionogram_{index:02d}.npz")
            truth = observed[index]
            scores.append(score(truth, modeled))
            nose_errors.append([abs(truth.nose(mode) - modeled.nose(mode))
                                for mode in (1, -1)])
        candidates.append({"name": name,
                           "ionograms": str(directory.relative_to(ROOT)),
                           "mean_score": float(np.mean([row["total"] for row in scores])),
                           "mean_return_distance": float(np.mean(
                               [row["return_distance"] for row in scores])),
                           "mean_ridge_score": float(np.mean([row["ridge"] for row in scores])),
                           "nose_mae_mhz_by_mode": np.mean(nose_errors, axis=0).tolist(),
                           "scores_by_profile": scores})
    chosen = min(candidates, key=lambda row: row["mean_score"])
    result = {"selection_uses_truth_density": False,
              "observed_ionograms": str(truth_dir.relative_to(ROOT)),
              "starting_retrieval": str(baseline_dir.relative_to(ROOT)),
              "stage": stage, "profile_indices": sorted(set(indices)),
              "selection_rule": "lowest mean 20-profile O/X ionogram score" if stage == "final"
              else "exploratory subset screen using stage-matched raw ionograms",
              "chosen_candidate": chosen["name"], "candidates": candidates}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"chosen_candidate": chosen["name"],
                      "scores": [{key: row[key] for key in
                                  ("name", "mean_score", "nose_mae_mhz_by_mode")}
                                 for row in candidates]}, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("raw", "final"), required=True)
    parser.add_argument("--candidates", nargs="+", required=True)
    parser.add_argument("--indices", nargs="+", type=int, default=list(range(1, 21)))
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    select(arguments.stage, arguments.candidates, arguments.indices, arguments.output)
