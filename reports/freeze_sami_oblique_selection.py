"""Freeze an oblique-ionogram choice before opening the full SAMI density field."""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

from evaluate_general_field_sami import profiles
from freeze_sami_monotone_selection import sha256

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "reports/data/sami3_monotone_oblique"
CASE = ROOT / "reports/data/sami3_wave_20170111_0600/oblique_wave_600km"
VERTICAL = ROOT / "reports/data/sami3_monotone_topside"
INDICES = [1, 5, 9, 12, 14, 17, 20]
INSITU_TOLERANCE_PERCENT = 0.5


def grid_for(name: str) -> Path:
    if name == "old_baseline":
        return CASE.parent / "pyiri_prior_grid.nc"
    if name == "old_fit":
        return CASE / "gain1_sigma60_height0.15_spline/grid.nc"
    if name == "varying_g04":
        return VERTICAL / "varying_g04/grid.nc"
    if name == "extended_anchor":
        return RUN / "extended_anchor/grid.nc"
    raise ValueError(f"Unknown candidate {name}")


def maximum_insitu_error_percent(grid: Path, samples: dict) -> float:
    altitude, modeled = profiles(grid, np.asarray(samples["lat_lon_deg"]))
    at_spacecraft = np.array([np.interp(samples["altitude_km"], altitude, row)
                              for row in modeled])
    return float(np.max(abs(100 * (at_spacecraft / np.asarray(samples["density_cm3"]) - 1))))


def main() -> None:
    os.umask(0o077)
    scores_path = RUN / "full_scores.json"
    scores = json.loads(scores_path.read_text())
    if scores["profile_indices"] != INDICES or scores["selection_uses_truth_density"]:
        raise ValueError("Expected seven ionogram-only scores")
    samples = json.loads((RUN / "insitu_1km.json").read_text())
    errors = {row["name"]: maximum_insitu_error_percent(grid_for(row["name"]), samples)
              for row in scores["candidates"]}
    feasible = [row for row in scores["candidates"]
                if errors[row["name"]] <= INSITU_TOLERANCE_PERCENT]
    if not feasible:
        raise ValueError("No candidate matches the 0.5% in-situ observations")
    name = min(feasible, key=lambda row: row["mean_total_equivalent_km"])["name"]
    grid = grid_for(name)
    row = next(item for item in scores["candidates"] if item["name"] == name)
    result = {
        "selection_uses_truth_density": False,
        "selection_uses_synthetic_insitu_density_samples": True,
        "selection_uses_full_3d_sami_density": False,
        "profile_indices": INDICES,
        "chosen_candidate": name,
        "chosen_mean_ionogram_score_km": row["mean_total_equivalent_km"],
        "grid": str(grid.relative_to(ROOT)),
        "grid_sha256": sha256(grid),
        "score_file_sha256": sha256(scores_path),
        "plan_file_sha256": sha256(RUN / "plan.json"),
        "insitu_max_absolute_error_percent_by_candidate": errors,
        "insitu_hard_tolerance_percent": INSITU_TOLERANCE_PERCENT,
        "choice_rule": "lowest mean uncapped O/X ionogram score among candidates matching all 1 km track samples within 0.5%; Doppler is a separately reported diagnostic",
    }
    output = RUN / "selection.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    output.chmod(0o600)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
