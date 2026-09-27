"""Refine an IRI D-ionogram fit with a local F2 peak correction.

The search and selection use observed O/X returns and the PyIRI background.
The separately generated IRI truth density is not opened here.

    python3 reports/fit_d_ionogram_peak.py start BASE_STATE WORKDIR BACKGROUND
    python3 reports/fit_d_ionogram_peak.py advance WORKDIR RESULTS_DIR
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from analyze_d_inverse_fit import _fitted_profile
from fit_d_ionogram import Ionogram
from fit_d_ionogram_topside import count_aware_score, _select_with_return_coverage


def _select_with_resolved_noses(state: dict) -> dict:
    """Require observed O/X cutoffs and return coverage before score ranking."""
    observed = Ionogram.read(Path(state["observed"]))
    counts = {mode: int(np.count_nonzero(observed.records[:, 1] == mode))
              for mode in (1, -1)}
    eligible = []
    for row in state["evaluations"]:
        predicted = Ionogram.read(Path(row["path"]))
        if all(predicted.nose(mode) is not None
               and abs(predicted.nose(mode) - observed.nose(mode)) < 0.05
               and abs(int(np.count_nonzero(predicted.records[:, 1] == mode)) - counts[mode])
               <= 0.10 * counts[mode] for mode in (1, -1)):
            eligible.append(row)
    if not eligible:
        raise RuntimeError("No candidate matched both resolved noses and return coverage")
    selected = min(eligible, key=lambda row: row["score"]["total"])
    state["selected_path"] = selected["path"]
    state["selection_rule"] = ("Each O/X nose within 0.05 MHz and each mode's "
                               "return count within 10% of observed; then minimum count-aware score")
    state["admissible_candidates"] = len(eligible)
    return selected


def start(base_state_path: Path, workdir: Path, background_path: Path) -> None:
    if (workdir / "state.json").exists():
        raise FileExistsError(f"Existing fit in {workdir}")
    base_state = json.loads(base_state_path.read_text())
    baseline = next(row for row in base_state["evaluations"]
                    if row["path"] == base_state["selected_path"])
    observed = Ionogram.read(Path(base_state["observed"]))
    o_nose = observed.nose(1)
    if o_nose is None or o_nose >= 10:
        raise ValueError("A resolved O-mode nose below 10 MHz is required")
    with np.load(background_path, allow_pickle=False) as background:
        latitude = float(background["tx_lat_deg"])
        longitude = float(background["tx_lon_deg"])
    _, profile = _fitted_profile(background_path, latitude, longitude, baseline)
    baseline_fof2 = float(0.00898 * np.sqrt(np.max(profile)))
    center = (o_nose / baseline_fof2) ** 2 - 1.0
    if not 0.005 <= center <= 0.1:
        raise ValueError(f"O-mode nose implies an unsuitable peak correction: {center:.3f}")
    # Bracket the observable O-mode nose with a deterministic one-dimensional
    # sweep. X mode and all ranges remain part of the selection score.
    factors = (0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0)
    candidates = []
    for task, factor in enumerate(factors, 1):
        candidates.append({"task": task,
                           **{key: baseline[key] for key in
                              ("density_scale", "hmf2_shift_km", "f2_width_scale",
                               "topside_width_ratio")},
                           "peak_perturbation_fraction": float(center * factor)})
    workdir.mkdir(parents=True, exist_ok=True, mode=0o700)
    state = {"observed": base_state["observed"],
             "base_state": str(base_state_path.resolve()),
             "background": str(background_path.resolve()),
             "search_model": "local_f2_peak_40km",
             "o_nose_mhz": o_nose,
             "baseline_fof2_mhz": baseline_fof2,
             "nose_center_fraction": center,
             "evaluations": base_state["evaluations"],
             "next_round": 1}
    (workdir / "state.json").write_text(json.dumps(state, indent=2) + "\n")
    (workdir / "candidates_round_1.json").write_text(json.dumps(candidates, indent=2) + "\n")
    print(json.dumps({"o_nose_mhz": o_nose, "baseline_fof2_mhz": baseline_fof2,
                      "nose_center_fraction": center, "candidate_count": len(candidates)}, indent=2))


def advance(workdir: Path, results_dir: Path) -> None:
    state_path = workdir / "state.json"
    state = json.loads(state_path.read_text())
    round_number = state["next_round"]
    if round_number not in (1, 2):
        raise ValueError("Peak refinement is complete")
    candidates = json.loads((workdir / f"candidates_round_{round_number}.json").read_text())
    observed = Ionogram.read(Path(state["observed"]))
    for candidate in candidates:
        path = results_dir / f"ionogram_{candidate['task']:02d}.npz"
        predicted = Ionogram.read(path)
        with np.load(path, allow_pickle=False) as data:
            actual = float(data["peak_perturbation_fraction"])
            actual_width = float(data["peak_width_km"])
        if abs(actual - candidate["peak_perturbation_fraction"]) > 1e-8:
            raise ValueError(f"Peak parameter does not match {path}")
        if abs(actual_width - candidate.get("peak_width_km", 40.0)) > 1e-8:
            raise ValueError(f"Peak width does not match {path}")
        if any(abs(getattr(predicted, key) - candidate[key]) > 1e-8 for key in
               ("density_scale", "hmf2_shift_km", "f2_width_scale", "topside_width_ratio")):
            raise ValueError(f"Base parameters do not match {path}")
        state["evaluations"].append({"path": str(path.resolve()),
                                     **{key: value for key, value in candidate.items()
                                        if key != "task"},
                                     "score": count_aware_score(observed, predicted)})
    selected = (_select_with_resolved_noses(state) if round_number == 2
                else _select_with_return_coverage(state))
    state["next_round"] = round_number + 1
    state_path.write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps({"candidate_count": len(state["evaluations"]),
                      "new_best": selected["path"], "score": selected["score"],
                      "admissible_candidates": state["admissible_candidates"]}, indent=2))


def propose_width(workdir: Path) -> None:
    state = json.loads((workdir / "state.json").read_text())
    if state["next_round"] != 2 or len(state["evaluations"]) != 248:
        raise ValueError("Width sweep requires a complete first peak sweep")
    path = workdir / "candidates_round_2.json"
    if path.exists():
        raise FileExistsError(path)
    original = json.loads(Path(state["base_state"]).read_text())
    baseline = next(row for row in original["evaluations"]
                    if row["path"] == original["selected_path"])
    center = state["nose_center_fraction"]
    candidates = []
    for width in (20.0, 30.0):
        for factor in (1.0, 1.25, 1.5, 1.75):
            candidates.append({"task": len(candidates) + 1,
                               **{key: baseline[key] for key in
                                  ("density_scale", "hmf2_shift_km", "f2_width_scale",
                                   "topside_width_ratio")},
                               "peak_perturbation_fraction": center * factor,
                               "peak_width_km": width})
    path.write_text(json.dumps(candidates, indent=2) + "\n")
    print(json.dumps({"candidate_count": len(candidates), "widths_km": [20, 30]}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    first = sub.add_parser("start")
    first.add_argument("base_state", type=Path)
    first.add_argument("workdir", type=Path)
    first.add_argument("background", type=Path)
    next_step = sub.add_parser("advance")
    next_step.add_argument("workdir", type=Path)
    next_step.add_argument("results_dir", type=Path)
    sub.add_parser("propose-width").add_argument("workdir", type=Path)
    args = parser.parse_args()
    if args.action == "start":
        start(args.base_state, args.workdir, args.background)
    elif args.action == "propose-width":
        propose_width(args.workdir)
    else:
        advance(args.workdir, args.results_dir)


if __name__ == "__main__":
    main()
