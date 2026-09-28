"""Select between the old and refitted wave candidates using ionograms only."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from fit_d_ionogram import Ionogram, score

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reports/data"
TRUTH = DATA / "lat_wave_continued_truth_ionograms"
CANDIDATES = {
    "previous_fit": DATA / "lat_wave_continued_previous_fit_ionograms",
    "new_fit": DATA / "lat_wave_continued_retrieved_ionograms",
}
OUTPUT = DATA / "lat_wave_continued_candidate_selection.json"


def select(truth_dir: Path = TRUTH, candidates: dict[str, Path] = CANDIDATES,
           output: Path = OUTPUT) -> None:
    scores = {name: [] for name in candidates}
    for index in range(1, 21):
        truth = Ionogram.read(truth_dir / f"ionogram_{index:02d}.npz")
        for name, directory in candidates.items():
            modeled = Ionogram.read(directory / f"ionogram_{index:02d}.npz")
            scores[name].append(score(truth, modeled))
    means = {name: float(np.mean([item["total"] for item in rows]))
             for name, rows in scores.items()}
    chosen = min(means, key=means.get)
    result = {
        "selection_uses_truth_density": False,
        "truth_ionograms": str(truth_dir.resolve().relative_to(ROOT)),
        "candidate_ionograms": {name: str(path.resolve().relative_to(ROOT))
                                 for name, path in candidates.items()},
        "selection_rule": "lowest mean 20-profile O/X ionogram score",
        "mean_ionogram_score": means,
        "score_by_profile": scores,
        "chosen_candidate": chosen,
    }
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"chosen_candidate": chosen, "mean_ionogram_score": means},
                     indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", type=Path, default=TRUTH)
    parser.add_argument("--previous", type=Path, default=CANDIDATES["previous_fit"])
    parser.add_argument("--new", type=Path, default=CANDIDATES["new_fit"])
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    select(args.truth, {"previous_fit": args.previous, "new_fit": args.new}, args.output)
