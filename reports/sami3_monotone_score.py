"""Score all accepted O/X returns without hiding large group-range errors."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from ionogram_metrics import Ionogram  # noqa: E402

CASE = ROOT / "reports/data/sami3_wave_20170111_0600"
OLD = ROOT / "reports/data/general_field_sami3_wave"
RUN = ROOT / "reports/data/sami3_monotone_topside"


def score_pair(observed_path: Path, modeled_path: Path) -> dict:
    observed = Ionogram.read(observed_path)
    modeled = Ionogram.read(modeled_path)
    if observed.settings != modeled.settings:
        raise ValueError("Candidate and truth generator settings differ")
    mode_results = {}
    for mode, label in ((1, "O"), (-1, "X")):
        a = observed.records[observed.records[:, 1] == mode]
        b = modeled.records[modeled.records[:, 1] == mode]
        bins = sorted(set(a[:, 0].astype(int)) | set(b[:, 0].astype(int)))
        errors = []
        missing = 0
        for index in bins:
            actual = a[a[:, 0] == index, 2]
            predicted = b[b[:, 0] == index, 2]
            if not len(actual) or not len(predicted):
                errors.append(150.0)
                missing += 1
            else:
                distance = abs(actual[:, None] - predicted[None, :])
                errors.append(float((np.mean(np.min(distance, axis=0))
                                     + np.mean(np.min(distance, axis=1))) / 2.0))
        nose_a, nose_b = observed.nose(mode), modeled.nose(mode)
        nose_error = (abs(nose_a - nose_b) if nose_a is not None and nose_b is not None
                      else 1.0)
        mode_results[label] = {
            "group_range_error_km": float(np.mean(errors)) if errors else 0.0,
            "missing_frequency_fraction": missing / len(bins) if bins else 0.0,
            "nose_error_mhz": float(nose_error),
            "frequencies_with_any_return": len(bins),
        }
    range_error = float(np.mean([item["group_range_error_km"]
                                 for item in mode_results.values()]))
    nose_error = float(np.mean([item["nose_error_mhz"]
                                for item in mode_results.values()]))
    return {"total_equivalent_km": range_error + 100.0 * nose_error,
            "mean_group_range_error_km": range_error,
            "mean_nose_error_mhz": nose_error,
            "by_mode": mode_results}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indices", nargs="+", type=int, required=True)
    parser.add_argument("--candidates", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    os.umask(0o077)
    rows = []
    for name in args.candidates:
        directory = (OLD / "forward/nose_init/ionograms" if name == "nose_init_previous" else
                     RUN / "forward" / name / "ionograms")
        parts = []
        for index in args.indices:
            parts.append({"index": index, **score_pair(
                CASE / "vertical_truth" / f"ionogram_{index:02d}.npz",
                directory / f"ionogram_{index:02d}.npz")})
        rows.append({"name": name,
                     "mean_total_equivalent_km": float(np.mean([
                         part["total_equivalent_km"] for part in parts])),
                     "mean_group_range_error_km": float(np.mean([
                         part["mean_group_range_error_km"] for part in parts])),
                     "by_profile": parts})
    result = {"selection_uses_truth_density": False,
              "profile_indices": args.indices,
              "score_definition": "mean symmetric nearest group-range error across every accepted O/X return at each occupied frequency, 150 km missing-frequency cost, plus 100 km per MHz mean O/X nose error; no range-error clipping",
              "candidates": rows,
              "winner": min(rows, key=lambda row: row["mean_total_equivalent_km"])["name"]}
    args.output.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.output.chmod(0o600)
    print(json.dumps({"winner": result["winner"],
                      "scores_km": {row["name"]: row["mean_total_equivalent_km"]
                                    for row in rows}}, indent=2))


if __name__ == "__main__":
    main()
