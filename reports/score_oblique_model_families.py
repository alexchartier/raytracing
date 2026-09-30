"""Select oblique shape candidates from O/X, Doppler, and in-situ observations."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from evaluate_general_field_sami import metrics, profiles  # noqa: E402
from ionogram_metrics import Ionogram  # noqa: E402
from oblique_wave_pass import doppler_score  # noqa: E402
from sami3_monotone_oblique_score import verify_geometry  # noqa: E402
from sami3_monotone_score import score_pair  # noqa: E402

RUN = ROOT / "reports/data/oblique_model_families"
CASES = ("nequick", "chapman", "iri2016")
DOPPLER_WEIGHT_KM = 20.0
INSITU_TOLERANCE_PERCENT = 0.5


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_for(case: str, name: str) -> Path:
    root = RUN / case
    return root / "forward" / name / "ionogram_01.npz"


def grid_for(case: str, name: str) -> Path:
    return RUN / case / name / "grid.nc"


def insitu_error(case: str, grid: Path) -> float:
    samples = json.loads((RUN / case / "insitu_1km.json").read_text())
    altitude, field = profiles(grid, np.asarray(samples["lat_lon_deg"]))
    prediction = np.array([np.interp(samples["altitude_km"], altitude, row)
                           for row in field])
    return float(np.max(abs(100 * (prediction /
                                 np.asarray(samples["density_cm3"]) - 1))))


def low_frequency_range_error(truth_path: Path, modeled_path: Path,
                              low_mhz: float = 2.5,
                              high_mhz: float = 4.5) -> float:
    """Symmetric accepted-return range error in a band, by occupied O/X mode."""
    observed = Ionogram.read(truth_path)
    modeled = Ionogram.read(modeled_path)
    if not np.array_equal(observed.frequencies, modeled.frequencies):
        raise ValueError("Different frequency grids")
    chosen_bins = set(np.flatnonzero((observed.frequencies >= low_mhz) &
                                    (observed.frequencies <= high_mhz)))
    per_mode = []
    for mode in (1, -1):
        actual = observed.records[observed.records[:, 1] == mode]
        predicted = modeled.records[modeled.records[:, 1] == mode]
        bins = sorted((set(actual[:, 0].astype(int)) |
                       set(predicted[:, 0].astype(int))) & chosen_bins)
        if not bins:
            continue
        errors = []
        for index in bins:
            a = actual[actual[:, 0] == index, 2]
            b = predicted[predicted[:, 0] == index, 2]
            if not len(a) or not len(b):
                errors.append(150.0)
            else:
                distance = abs(a[:, None] - b[None, :])
                errors.append(float((np.mean(np.min(distance, axis=0)) +
                                     np.mean(np.min(distance, axis=1))) / 2))
        per_mode.append(float(np.mean(errors)))
    if not per_mode:
        raise ValueError("No accepted returns in the frequency band")
    return float(np.mean(per_mode))


def score_case(case: str, candidates: list[str]) -> None:
    if case not in CASES:
        raise ValueError(case)
    truth = path_for(case, "truth")
    rows = []
    for name in candidates:
        modeled = path_for(case, name)
        verify_geometry(truth, modeled, 1)
        ionogram = score_pair(truth, modeled)
        doppler = doppler_score(truth, modeled)
        error = insitu_error(case, grid_for(case, name))
        rows.append({"name": name,
                     "ionogram_equivalent_km": ionogram["total_equivalent_km"],
                     "group_range_error_km": ionogram["mean_group_range_error_km"],
                     "low_band_2p5_to_4p5_error_km": low_frequency_range_error(
                         truth, modeled),
                     "nose_error_mhz": ionogram["mean_nose_error_mhz"],
                     "by_mode": ionogram["by_mode"],
                     "doppler_1hz_mismatch": doppler,
                     "insitu_max_absolute_error_percent": error,
                     "combined_equivalent_km": ionogram["total_equivalent_km"]
                     + DOPPLER_WEIGHT_KM * doppler})
    feasible = [row for row in rows
                if row["insitu_max_absolute_error_percent"] <= INSITU_TOLERANCE_PERCENT]
    if not feasible:
        raise ValueError("No candidate satisfies the 800 km observation")
    chosen = min(feasible, key=lambda row: row["combined_equivalent_km"])
    result = {
        "case": case,
        "uses_full_truth_density": False,
        "uses_synthetic_800km_track_density": True,
        "ionogram_truth": str(truth.relative_to(ROOT)),
        "score": "uncapped symmetric nearest O/X group range, 150 km missing-frequency cost, 100 km/MHz mean nose cost; plus 20 km times 1 Hz Doppler mismatch",
        "doppler_weight_equivalent_km": DOPPLER_WEIGHT_KM,
        "insitu_tolerance_percent": INSITU_TOLERANCE_PERCENT,
        "candidates": rows,
        "lowest_combined_candidate": chosen["name"],
    }
    dest = RUN / case / "scores.json"
    dest.write_text(json.dumps(result, indent=2) + "\n")
    dest.chmod(0o600)
    print(json.dumps({"case": case,
                      "scores": {row["name"]: round(row["combined_equivalent_km"], 2)
                                 for row in rows},
                      "lowest": chosen["name"]}, indent=2))


def freeze(case: str) -> None:
    root = RUN / case
    scores_path = root / "scores.json"
    scores = json.loads(scores_path.read_text())
    if scores["case"] != case or scores["uses_full_truth_density"]:
        raise ValueError("Scores were not made from observables only")
    name = scores["lowest_combined_candidate"]
    grid = grid_for(case, name)
    result = {
        "case": case, "chosen_candidate": name,
        "uses_full_truth_density": False,
        "uses_synthetic_800km_track_density": True,
        "scores_sha256": digest(scores_path),
        "selected_grid_sha256": digest(grid),
        "selected_grid": str(grid.relative_to(ROOT)),
        "selection_rule": "minimum combined O/X plus 1 Hz Doppler score among fields within 0.5% of all 800 km observations",
    }
    path = root / "selection.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    path.chmod(0o600)
    print(json.dumps(result, indent=2))


def evaluate(case: str) -> None:
    root = RUN / case
    scores_path = root / "scores.json"
    scores = json.loads(scores_path.read_text())
    selection = json.loads((root / "selection.json").read_text())
    name = selection["chosen_candidate"]
    grid = grid_for(case, name)
    if (selection["uses_full_truth_density"] or
            selection["scores_sha256"] != digest(scores_path) or
            selection["selected_grid_sha256"] != digest(grid) or
            scores["lowest_combined_candidate"] != name):
        raise ValueError("Candidate choice must be frozen before density evaluation")
    manifest = json.loads((root / "manifest.json").read_text())
    source = manifest["profiles"][0]
    with np.load(path_for(case, "truth"), allow_pickle=False) as ray:
        midpoint = (float(ray["tx_lat_deg"]) + float(ray["rx_lat_deg"])) / 2
    locations = np.array([[source["latitude_deg"], source["longitude_deg"]],
                          [midpoint, source["longitude_deg"]]])
    altitude, truth = profiles(root / "truth_grid.nc", locations)
    _, start = profiles(grid_for(case, "start"), locations)
    _, selected = profiles(grid, locations)
    result = {
        "case": case, "truth_model": manifest["density_model"],
        "selected": name, "truth_density_read_after_selection": True,
        "positions": ["transmitter", "link_midpoint"],
        "starting_density_metrics": metrics(altitude, truth, start),
        "selected_density_metrics": metrics(altitude, truth, selected),
        "starting_ionogram_score": next(row for row in scores["candidates"]
                                         if row["name"] == "start"),
        "selected_ionogram_score": next(row for row in scores["candidates"]
                                         if row["name"] == name),
    }
    dest = root / "evaluation.json"
    dest.write_text(json.dumps(result, indent=2) + "\n")
    dest.chmod(0o600)
    print(json.dumps({"case": case, "selected": name,
                      "density": result["selected_density_metrics"]}, indent=2))


def diagnose(case: str) -> None:
    """Record low-frequency detail after the observable selection is frozen."""
    root = RUN / case
    selection = json.loads((root / "selection.json").read_text())
    scores_path = root / "scores.json"
    if (selection["uses_full_truth_density"] or
            selection["scores_sha256"] != digest(scores_path)):
        raise ValueError("Freeze the observable selection first")
    truth = path_for(case, "truth")
    bands = ((2.0, 2.4), (2.0, 4.5), (2.5, 4.5))

    def band_result(name: str, low: float, high: float) -> float | None:
        try:
            return low_frequency_range_error(truth, path_for(case, name),
                                             low, high)
        except ValueError as error:
            if str(error) != "No accepted returns in the frequency band":
                raise
            return None

    result = {"case": case, "selection_frozen": True,
              "uses_full_truth_density": False,
              "method": "mean symmetric nearest O/X group range, 150 km for a missing frequency; only modes with occupied bins contribute",
              "bands": {f"{low:.1f}_to_{high:.1f}_mhz": {
                  name: band_result(name, low, high)
                  for name in ("start", selection["chosen_candidate"])}
                  for low, high in bands}}
    target = root / "low_frequency_diagnostic.json"
    target.write_text(json.dumps(result, indent=2) + "\n")
    target.chmod(0o600)
    print(json.dumps(result, indent=2))


def main() -> None:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("score", "freeze", "evaluate", "diagnose"))
    parser.add_argument("--case", required=True, choices=CASES)
    parser.add_argument("--candidates", nargs="+")
    args = parser.parse_args()
    if args.action == "score":
        if not args.candidates:
            parser.error("--candidates required when scoring")
        score_case(args.case, args.candidates)
    elif args.action == "freeze":
        freeze(args.case)
    elif args.action == "diagnose":
        diagnose(args.case)
    else:
        evaluate(args.case)


if __name__ == "__main__":
    main()
