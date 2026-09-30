"""Freeze the seven-ionogram candidate choice before reading SAMI density truth."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "reports/data/sami3_monotone_topside"
INDICES = [1, 5, 9, 12, 14, 17, 20]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    os.umask(0o077)
    score_path = RUN / "full_scores.json"
    plan_path = RUN / "plan.json"
    scores = json.loads(score_path.read_text())
    plan = json.loads(plan_path.read_text())
    if scores["profile_indices"] != INDICES or scores["selection_uses_truth_density"]:
        raise ValueError("Expected seven ionogram-only scoring profiles")
    names = {row["name"] for row in scores["candidates"]}
    if names != {"nose_init_previous", "varying_g04", "varying_g07"}:
        raise ValueError("Expected baseline and the two pilot finalists")
    winner = min(scores["candidates"],
                 key=lambda row: row["mean_total_equivalent_km"])
    if winner["name"] == "nose_init_previous":
        raise ValueError("No new candidate improves the previous field")
    grid = ROOT / next(row["grid"] for row in plan["candidates"]
                       if row["name"] == winner["name"])
    result = {
        "chosen_candidate": winner["name"],
        "chosen_mean_total_equivalent_km": winner["mean_total_equivalent_km"],
        "profile_indices": INDICES,
        "selection_uses_truth_density": False,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "score_file_sha256": sha256(score_path),
        "plan_file_sha256": sha256(plan_path),
        "selected_grid_sha256": sha256(grid),
    }
    output = RUN / "selection.json"
    if output.exists():
        raise FileExistsError("Selection is already frozen")
    output.write_text(json.dumps(result, indent=2) + "\n")
    output.chmod(0o600)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
