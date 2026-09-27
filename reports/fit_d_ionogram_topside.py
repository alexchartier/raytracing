"""Extend a completed IRI D-ionogram fit with a separate topside width.

The withheld truth density is never read. Candidates and selection use only
accepted ionogram returns.

    python3 reports/fit_d_ionogram_topside.py start BASE_STATE WORKDIR
    python3 reports/fit_d_ionogram_topside.py advance WORKDIR RESULTS_DIR

Each advance ingests a 40-member population. Two populations are proposed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from fit_d_ionogram import Ionogram, score

BOUNDS = np.array([[0.25, 1.40], [-100.0, 100.0], [0.60, 1.85], [0.55, 1.35]])
POPULATION = 40


def count_aware_score(observed: Ionogram, predicted: Ionogram) -> dict[str, float]:
    """Fit range, nose, ridge, and mode-wise accepted-return counts."""
    base = score(observed, predicted)
    count_errors = []
    for mode in (1, -1):
        observed_count = int(np.count_nonzero(observed.records[:, 1] == mode))
        predicted_count = int(np.count_nonzero(predicted.records[:, 1] == mode))
        count_errors.append(min(abs(observed_count - predicted_count)
                                / max(observed_count, 1), 1.0))
    count_penalty = float(np.mean(count_errors))
    return {**base, "legacy_total": base["total"],
            "count_penalty": count_penalty,
            "total": float(.50 * base["return_distance"] + .15 * base["nose"]
                           + .20 * base["ridge"] + .15 * count_penalty)}


def _rescore_history(state: dict) -> None:
    observed = Ionogram.read(Path(state["observed"]))
    for row in state["evaluations"]:
        row["score"] = count_aware_score(observed, Ionogram.read(Path(row["path"])))


def _select_with_return_coverage(state: dict) -> dict:
    """Require each mode's return count within 10% before score ranking."""
    observed = Ionogram.read(Path(state["observed"]))
    truth_counts = {mode: int(np.count_nonzero(observed.records[:, 1] == mode))
                    for mode in (1, -1)}
    eligible = []
    for row in state["evaluations"]:
        predicted = Ionogram.read(Path(row["path"]))
        counts = {mode: int(np.count_nonzero(predicted.records[:, 1] == mode))
                  for mode in (1, -1)}
        if all(abs(counts[mode] - truth_counts[mode]) <= .10 * truth_counts[mode]
               for mode in (1, -1)):
            eligible.append(row)
    if not eligible:
        raise RuntimeError("No candidate met the per-mode return coverage rule")
    selected = min(eligible, key=lambda row: row["score"]["total"])
    state["selected_path"] = selected["path"]
    state["selection_rule"] = "Each O/X return count within 10% of observed; then minimum count-aware score"
    state["admissible_candidates"] = len(eligible)
    return selected


def _point(row: dict) -> np.ndarray:
    return np.array([row["density_scale"], row["hmf2_shift_km"],
                     row.get("f2_width_scale", 1.0), row.get("topside_width_ratio", 1.0)])


def _candidate(task: int, point: np.ndarray) -> dict:
    return {"task": task, "density_scale": float(point[0]),
            "hmf2_shift_km": float(point[1]), "f2_width_scale": float(point[2]),
            "topside_width_ratio": float(point[3])}


def _unique(points: list[np.ndarray], history: list[dict], rng: np.random.Generator) -> list[dict]:
    previous = np.array([_point(row) for row in history])
    result = []
    tolerance = np.array([1e-5, 1e-3, 1e-4, 1e-4])
    for proposed in points:
        point = np.clip(proposed, BOUNDS[:, 0], BOUNDS[:, 1])
        for _ in range(200):
            duplicate = np.any(np.all(np.abs(previous - point) < tolerance, axis=1))
            duplicate |= any(np.all(np.abs(_point(row) - point) < tolerance) for row in result)
            if not duplicate:
                break
            point = np.clip(point + rng.normal(size=4) * [.003, .5, .01, .01],
                            BOUNDS[:, 0], BOUNDS[:, 1])
        else:
            raise RuntimeError("Could not propose unique topside candidates")
        result.append(_candidate(len(result) + 1, point))
    return result


def _save(workdir: Path, state: dict, candidates: list[dict]) -> None:
    workdir.mkdir(parents=True, exist_ok=True, mode=0o700)
    (workdir / "state.json").write_text(json.dumps(state, indent=2) + "\n")
    (workdir / f"candidates_round_{state['next_round']}.json").write_text(
        json.dumps(candidates, indent=2) + "\n")
    print(json.dumps({"next_round": state["next_round"], "candidate_count": len(candidates),
                      "best_score": min(row["score"]["total"] for row in state["evaluations"])},
                     indent=2))


def _propose(state: dict) -> list[dict]:
    evaluations = state["evaluations"]
    ranked = sorted(evaluations, key=lambda row: row["score"]["total"])
    best = _point(ranked[0])
    round_number = state["next_round"]
    rng = np.random.default_rng(20260928 + round_number)
    points = []
    if round_number == 1:
        # Cover the density/height/width neighborhood while testing narrower
        # and wider topsides. The old best, with ratio 1, remains in the pool.
        sample = qmc.LatinHypercube(d=4, seed=20260928).random(28)
        for u in sample:
            points.append(np.array([.39 + .14*u[0], -83 + 39*u[1],
                                    1.17 + .58*u[2], .63 + .58*u[3]]))
        for _ in range(12):
            points.append(best + rng.normal(size=4) * [.022, 6.0, .09, .16])
    elif round_number == 2:
        elite = np.array([_point(row) for row in ranked[:8]])
        for i in range(32):
            center = elite[i % 8] if i % 4 == 0 else best
            points.append(center + rng.normal(size=4) * [.016, 4.5, .065, .07])
        for _ in range(8):
            points.append(best + rng.normal(size=4) * [.034, 9.0, .14, .14])
    else:
        raise ValueError("Only two topside populations are configured")
    return _unique(points, evaluations, rng)


def start(base_state: Path, workdir: Path) -> None:
    if (workdir / "state.json").exists():
        raise FileExistsError(f"Existing fit in {workdir}")
    original = json.loads(base_state.read_text())
    if len(original["evaluations"]) != 160 or original["next_round"] != 4:
        raise ValueError("Expected a completed 160-candidate fit")
    state = {"observed": original["observed"], "base_state": str(base_state.resolve()),
             "search_model": "scale_height_width_topside", "next_round": 1,
             "evaluations": original["evaluations"]}
    _rescore_history(state)
    _save(workdir, state, _propose(state))


def reweight(workdir: Path) -> None:
    """Apply count-aware scoring before the second population is submitted."""
    state_path = workdir / "state.json"
    state = json.loads(state_path.read_text())
    if state["next_round"] != 2 or len(state["evaluations"]) != 200:
        raise ValueError("Reweighting requires the complete first topside population")
    _rescore_history(state)
    _save(workdir, state, _propose(state))


def advance(workdir: Path, results_dir: Path) -> None:
    state_path = workdir / "state.json"
    state = json.loads(state_path.read_text())
    round_number = state["next_round"]
    if round_number not in (1, 2):
        raise ValueError("Topside fit is complete")
    candidates = json.loads((workdir / f"candidates_round_{round_number}.json").read_text())
    if len(candidates) != POPULATION:
        raise ValueError("Expected 40 candidates")
    observed = Ionogram.read(Path(state["observed"]))
    for candidate in candidates:
        path = results_dir / f"ionogram_{candidate['task']:02d}.npz"
        ionogram = Ionogram.read(path)
        actual = np.array([ionogram.density_scale, ionogram.hmf2_shift_km,
                           ionogram.f2_width_scale, ionogram.topside_width_ratio])
        if not np.allclose(actual, _point(candidate), rtol=0, atol=1e-8):
            raise ValueError(f"Candidate parameters do not match {path}")
        state["evaluations"].append({"path": str(path.resolve()),
                                      **{key: candidate[key] for key in
                                         ("density_scale", "hmf2_shift_km",
                                          "f2_width_scale", "topside_width_ratio")},
                                      "score": count_aware_score(observed, ionogram)})
    state["next_round"] += 1
    if round_number == 2:
        selected = _select_with_return_coverage(state)
        state_path.write_text(json.dumps(state, indent=2) + "\n")
        best = min(state["evaluations"], key=lambda row: row["score"]["total"])
        print(json.dumps({"complete": True, "evaluations": len(state["evaluations"]),
                          "lowest_score": best, "selected": selected}, indent=2))
    else:
        _save(workdir, state, _propose(state))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    first = sub.add_parser("start")
    first.add_argument("base_state", type=Path)
    first.add_argument("workdir", type=Path)
    next_step = sub.add_parser("advance")
    next_step.add_argument("workdir", type=Path)
    next_step.add_argument("results_dir", type=Path)
    sub.add_parser("reweight").add_argument("workdir", type=Path)
    sub.add_parser("select").add_argument("workdir", type=Path)
    args = parser.parse_args()
    if args.action == "start":
        start(args.base_state, args.workdir)
    elif args.action == "reweight":
        reweight(args.workdir)
    elif args.action == "select":
        workdir = args.workdir
        path = workdir / "state.json"
        state = json.loads(path.read_text())
        if state["next_round"] != 3 or len(state["evaluations"]) != 240:
            raise ValueError("Selection requires all 240 evaluations")
        selected = _select_with_return_coverage(state)
        path.write_text(json.dumps(state, indent=2) + "\n")
        print(json.dumps({"selected": selected,
                          "admissible_candidates": state["admissible_candidates"]}, indent=2))
    else:
        advance(args.workdir, args.results_dir)


if __name__ == "__main__":
    main()
