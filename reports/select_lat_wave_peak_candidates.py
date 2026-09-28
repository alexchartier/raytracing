"""Rank 20-profile peak-correction candidates using O/X ionograms only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from fit_d_ionogram import Ionogram, score

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reports/data"
MANIFEST = DATA / "lat_wave_peak_candidates/candidates.json"


def rank(observed_dir: Path, baseline_dir: Path, candidate_stage: str,
         output: Path, candidate_names: list[str] | None = None) -> dict:
    configuration = json.loads(MANIFEST.read_text())
    known = {row["name"] for row in configuration["candidates"]}
    if candidate_names is not None and (not candidate_names or
                                        any(name not in known for name in candidate_names)):
        raise ValueError("Candidate names must come from the candidate manifest")
    candidates = {"baseline": baseline_dir}
    for row in configuration["candidates"]:
        if candidate_names is not None and row["name"] not in candidate_names:
            continue
        candidates[row["name"]] = DATA / "lat_wave_peak_candidates" / row["name"] / candidate_stage
    observed = [Ionogram.read(observed_dir / f"ionogram_{index:02d}.npz")
                for index in range(1, 21)]
    rows = []
    for name, directory in candidates.items():
        per_profile = []
        nose_mae = []
        for index, truth in enumerate(observed, 1):
            modeled = Ionogram.read(directory / f"ionogram_{index:02d}.npz")
            per_profile.append(score(truth, modeled))
            nose_mae.append([abs(truth.nose(mode) - modeled.nose(mode))
                             for mode in (1, -1)])
        rows.append({"name": name, "ionograms": str(directory.resolve().relative_to(ROOT)),
                     "mean_score": float(np.mean([item["total"] for item in per_profile])),
                     "nose_mae_mhz_by_mode": np.mean(nose_mae, axis=0).tolist(),
                     "scores_by_profile": per_profile})
    chosen = min(rows, key=lambda item: item["mean_score"])
    result = {
        "selection_uses_truth_density": False,
        "observed_ionograms": str(observed_dir.resolve().relative_to(ROOT)),
        "baseline_ionograms": str(baseline_dir.resolve().relative_to(ROOT)),
        "candidate_stage": candidate_stage,
        "selection_rule": "lowest mean 20-profile O/X ionogram score",
        "candidate_names_evaluated": list(candidates),
        "chosen_candidate": chosen["name"],
        "candidates": rows,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"chosen_candidate": chosen["name"],
                      "scores": [{key: row[key] for key in
                                  ("name", "mean_score", "nose_mae_mhz_by_mode")}
                                 for row in rows]}, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observed", type=Path, required=True)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--stage", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--candidates", nargs="+",
                        help="Candidate variants advanced to this stage")
    args = parser.parse_args()
    rank(args.observed, args.baseline, args.stage, args.output, args.candidates)


if __name__ == "__main__":
    main()
