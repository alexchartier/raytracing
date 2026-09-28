"""Rank low-Doppler peak candidates using observed O/X ionograms only."""

from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np

from fit_d_ionogram import Ionogram, score

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reports/data"
ROUND = DATA / "lat_wave_doppler_peak_round3"
LIMIT_HZ = 15.0


def gated_nose(path: Path, mode: int) -> float | None:
    with np.load(path, allow_pickle=False) as data:
        records = np.asarray(data["records"], dtype=float)
        frequencies = np.asarray(data["frequencies_mhz"], dtype=float)
        doppler = np.asarray(data["spacecraft_doppler_hz"], dtype=float)
    selected = records[(records[:, 1] == mode) & (abs(doppler) <= LIMIT_HZ)]
    return float(frequencies[int(np.max(selected[:, 0]))]) if len(selected) else None


def rank(stage: str, indices: list[int], names: list[str], output: Path) -> dict:
    known = {row["name"] for row in
             json.loads((ROUND / "candidates.json").read_text())["candidates"]}
    if not names or any(name not in known for name in names):
        raise ValueError("Choose candidate names from candidates.json")
    if stage == "final" and sorted(set(indices)) != list(range(1, 21)):
        raise ValueError("Final selection requires all 20 profiles")
    if not indices or any(index < 1 or index > 20 for index in indices):
        raise ValueError("Invalid profile indices")
    # The original raw truth files predate Doppler storage, so the exploratory
    # raw screen uses final observed rays against each candidate's raw rays.
    truth_dir = DATA / "lat_wave_final_truth_ionograms"
    directories = {"baseline": DATA / "lat_wave_joint_round2/peak_full" / stage}
    directories.update({name: ROUND / name / stage for name in names})
    rows = []
    for name, directory in directories.items():
        scores = []
        gate_abs_error = []
        for index in indices:
            filename = f"ionogram_{index:02d}.npz"
            truth = Ionogram.read(truth_dir / filename)
            modeled = Ionogram.read(directory / filename)
            scores.append(score(replace(truth, settings=modeled.settings)
                                if stage == "raw" else truth, modeled))
            gate_abs_error.append([abs(a - b) if a is not None and b is not None
                                   else .8 for mode in (1, -1)
                                   for a, b in [(gated_nose(truth_dir / filename, mode),
                                                 gated_nose(directory / filename, mode))]])
        standard = float(np.mean([row["total"] for row in scores]))
        gated = float(np.mean(np.minimum(np.asarray(gate_abs_error) / .8, 1.0)))
        rows.append({"name": name, "ionograms": str(directory.relative_to(ROOT)),
                     "standard_score": standard,
                     "gated_cutoff_score": gated,
                     "combined_score": .7 * standard + .3 * gated,
                     "gated_cutoff_mae_mhz_by_mode": np.mean(gate_abs_error, axis=0).tolist(),
                     "scores_by_profile": scores})
    chosen = min(rows, key=lambda row: row["combined_score"])
    result = {"selection_uses_truth_density": False,
              "observed_ionograms": str(truth_dir.relative_to(ROOT)),
              "stage": stage, "profile_indices": sorted(set(indices)),
              "raw_screen_uses_final_observations": stage == "raw",
              "doppler_limit_hz": LIMIT_HZ,
              "selection_rule": "lowest mean (0.7 standard O/X score + 0.3 low-Doppler cutoff score)",
              "chosen_candidate": chosen["name"], "candidates": rows}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"chosen_candidate": chosen["name"],
                      "scores": [{key: row[key] for key in
                                  ("name", "standard_score", "gated_cutoff_score",
                                   "combined_score", "gated_cutoff_mae_mhz_by_mode")}
                                 for row in rows]}, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("raw", "final"), required=True)
    parser.add_argument("--indices", nargs="+", type=int, default=list(range(1, 21)))
    parser.add_argument("--candidates", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rank(args.stage, args.indices, args.candidates, args.output)
