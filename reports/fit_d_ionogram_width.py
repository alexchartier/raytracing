"""Explore F2 profile width after a 40+40 scale/height ionogram search.

This is an exploratory model extension. The candidate score uses only accepted
ionogram returns; a truth density grid is never read by this program.

    python3 reports/fit_d_ionogram_width.py start FIT_DIR
    python3 reports/fit_d_ionogram_width.py advance FIT_DIR RESULTS_DIR

``start`` replaces the unrun round-2 candidate proposal with 40 three-parameter
proposals. Two ``advance`` calls score rounds 2 and 3, for 160 total ionograms.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import qmc

from fit_d_ionogram import Ionogram, score

BOUNDS = np.array([[0.25, 1.40], [-100.0, 100.0], [0.60, 1.85]])


def _parameters(row: dict) -> np.ndarray:
    return np.array([row["density_scale"], row["hmf2_shift_km"], row.get("f2_width_scale", 1.0)])


def _proposal(task: int, point: np.ndarray) -> dict:
    return {"task": task, "density_scale": float(point[0]),
            "hmf2_shift_km": float(point[1]), "f2_width_scale": float(point[2])}


def _unique(points: list[np.ndarray], previous: np.ndarray, rng: np.random.Generator) -> list[dict]:
    result = []
    for point in points:
        point = np.clip(point, BOUNDS[:, 0], BOUNDS[:, 1])
        for _ in range(200):
            if (not np.any(np.all(np.abs(previous - point) < [1e-5, 1e-3, 1e-4], axis=1))
                    and not any(np.all(np.abs(_parameters(row) - point) < [1e-5, 1e-3, 1e-4])
                                for row in result)):
                break
            point = np.clip(point + rng.normal(size=3) * [.003, .5, .01],
                            BOUNDS[:, 0], BOUNDS[:, 1])
        else:
            raise RuntimeError("Could not produce unique width candidates")
        result.append(_proposal(len(result) + 1, point))
    return result


def _save(state_path: Path, state: dict, candidates: list[dict], round_number: int) -> None:
    state_path.write_text(json.dumps(state, indent=2) + "\n")
    candidate_path = state_path.parent / f"candidates_round_{round_number}.json"
    candidate_path.write_text(json.dumps(candidates, indent=2) + "\n")
    best = min(state["evaluations"], key=lambda row: row["score"]["total"])
    print(json.dumps({"round": round_number, "candidate_count": len(candidates),
                      "best": best}, indent=2))


def start(workdir: Path) -> None:
    state_path = workdir / "state.json"
    state = json.loads(state_path.read_text())
    if state["next_round"] != 2 or not 70 <= len(state["evaluations"]) <= 80:
        raise ValueError("Width search starts after at least 30 first-round candidates")
    if state.get("search_model") == "scale_height_width":
        raise ValueError("Width search was already started")
    best = min(state["evaluations"], key=lambda row: row["score"]["total"])
    center = _parameters(best)
    rng = np.random.default_rng(20260926)
    points = []
    # A broad Latin hypercube prevents the 2D optimum from locking the
    # three-parameter search into one profile shape.
    sample = qmc.LatinHypercube(d=3, seed=20260926).random(32)
    for u, v, w in sample:
        points.append(np.array([center[0] + .26 * (u - .5),
                                center[1] + 70.0 * (v - .5),
                                .70 + 1.05 * w]))
    for _ in range(8):
        points.append(center + rng.normal(size=3) * [.08, 18.0, .25])
    previous = np.array([_parameters(row) for row in state["evaluations"]])
    candidates = _unique(points, previous, rng)
    state["search_model"] = "scale_height_width"
    _save(state_path, state, candidates, 2)


def adopt(workdir: Path, proposal_path: Path, seed_evaluations: int) -> None:
    """Attach an early, ionogram-only proposal to the completed 80-case state."""
    state_path = workdir / "state.json"
    state = json.loads(state_path.read_text())
    if state["next_round"] != 2 or len(state["evaluations"]) != 80:
        raise ValueError("Adoption requires a complete first refinement round")
    if not 70 <= seed_evaluations <= 80:
        raise ValueError("Expected at least 70 evaluations at proposal time")
    candidates = json.loads(proposal_path.read_text())
    if len(candidates) != 40 or [row["task"] for row in candidates] != list(range(1, 41)):
        raise ValueError("Expected 40 numbered width candidates")
    if any("f2_width_scale" not in row for row in candidates):
        raise ValueError("Missing width parameter")
    state["search_model"] = "scale_height_width"
    state["round_2_proposal_evaluations"] = seed_evaluations
    _save(state_path, state, candidates, 2)


def advance(workdir: Path, results_dir: Path) -> None:
    state_path = workdir / "state.json"
    state = json.loads(state_path.read_text())
    round_number = state["next_round"]
    if state.get("search_model") != "scale_height_width" or round_number not in (2, 3):
        raise ValueError("Expected an active width search in round 2 or 3")
    candidates = json.loads((workdir / f"candidates_round_{round_number}.json").read_text())
    observed = Ionogram.read(Path(state["observed"]))
    new = []
    for candidate in candidates:
        path = results_dir / f"ionogram_{candidate['task']:02d}.npz"
        ionogram = Ionogram.read(path)
        actual = np.array([ionogram.density_scale, ionogram.hmf2_shift_km, ionogram.f2_width_scale])
        if not np.allclose(actual, _parameters(candidate), rtol=0, atol=[1e-8, 1e-8, 1e-8]):
            raise ValueError(f"Candidate parameters do not match {path}")
        new.append({"path": str(path.resolve()), "density_scale": actual[0],
                    "hmf2_shift_km": actual[1], "f2_width_scale": actual[2],
                    "score": score(observed, ionogram)})
    state["evaluations"].extend(new)
    state["next_round"] += 1
    if round_number == 3:
        state_path.write_text(json.dumps(state, indent=2) + "\n")
        best = min(state["evaluations"], key=lambda row: row["score"]["total"])
        print(json.dumps({"complete": True, "evaluations": len(state["evaluations"]),
                          "best": best}, indent=2))
        return
    ranked = sorted(state["evaluations"], key=lambda row: row["score"]["total"])
    best = _parameters(ranked[0])
    elite = np.array([_parameters(row) for row in ranked[:8]])
    rng = np.random.default_rng(20260927)
    points = []
    for i in range(32):
        center = elite[i % len(elite)] if i % 4 == 0 else best
        points.append(center + rng.normal(size=3) * [.025, 6.0, .10])
    for _ in range(8):
        points.append(best + rng.normal(size=3) * [.065, 14.0, .22])
    previous = np.array([_parameters(row) for row in state["evaluations"]])
    _save(state_path, state, _unique(points, previous, rng), 3)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("start").add_argument("workdir", type=Path)
    attach = sub.add_parser("adopt")
    attach.add_argument("workdir", type=Path)
    attach.add_argument("proposal_path", type=Path)
    attach.add_argument("seed_evaluations", type=int)
    step = sub.add_parser("advance")
    step.add_argument("workdir", type=Path)
    step.add_argument("results_dir", type=Path)
    args = parser.parse_args()
    if args.action == "start":
        start(args.workdir)
    elif args.action == "adopt":
        adopt(args.workdir, args.proposal_path, args.seed_evaluations)
    else:
        advance(args.workdir, args.results_dir)


if __name__ == "__main__":
    main()
