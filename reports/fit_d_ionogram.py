"""Population fit of a D-generator vertical ionogram (density and F2 height).

The input is an observed accepted-return NPZ, never the generator's truth
parameters. ``init`` scores a library of 40 baseline ionograms and proposes
40 new candidates. ``advance`` scores each completed batch and proposes the
next one. After three advances, the four 40-member populations are complete.

    python3 reports/fit_d_ionogram.py init OBSERVED LIBRARY RUN_DIR
    python3 reports/fit_d_ionogram.py advance RUN_DIR RESULTS_DIR

Candidates are in RUN_DIR/candidates_round_N.json. Generate each candidate
with generate_synthetic_truth_returns.py using its density_scale and
hmf2_shift_km, and name outputs ionogram_01.npz ... ionogram_40.npz.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

BOUNDS = np.array([[0.85, 1.35], [-40.0, 40.0]])
POPULATION = 40
REFINEMENTS = 3


@dataclass(frozen=True)
class Ionogram:
    frequencies: np.ndarray
    records: np.ndarray
    density_scale: float | None
    hmf2_shift_km: float | None
    settings: tuple

    @classmethod
    def read(cls, path: Path) -> "Ionogram":
        with np.load(path, allow_pickle=False) as data:
            frequencies = np.asarray(data["frequencies_mhz"], dtype=float)
            records = np.asarray(data["records"], dtype=float)
            density = float(data["density_scale"]) if "density_scale" in data else None
            shift = float(data["hmf2_shift_km"]) if "hmf2_shift_km" in data else 0.0
            settings = (str(data["method"]), str(data["vertical_fan_layout"]),
                        float(data["vertical_outer_ray_fraction"]),
                        int(data["vertical_guard_seed_limit"]),
                        float(data["homing_tolerance_m"]))
        if records.ndim != 2 or records.shape[1] < 3:
            raise ValueError(f"Bad records in {path}")
        if len(frequencies) != 81 or not np.allclose(frequencies, np.arange(2, 10.0001, .1)):
            raise ValueError(f"Expected the complete 2–10 MHz, 100 kHz sweep in {path}")
        if len(records) and (np.min(records[:, 0]) < 0 or np.max(records[:, 0]) >= len(frequencies)):
            raise ValueError(f"Bad frequency indices in {path}")
        if settings != ("adaptive", "equal_area_guarded", .5, 4, 1000.0):
            raise ValueError(f"Expected the D generator and 1 km homing gate in {path}: {settings}")
        return cls(frequencies, records, density, shift, settings)

    def nose(self, mode: int) -> float | None:
        indices = self.records[self.records[:, 1] == mode, 0].astype(int)
        return float(self.frequencies[np.max(indices)]) if len(indices) else None

    def ridge(self, mode: int) -> dict[int, float]:
        records = self.records[self.records[:, 1] == mode]
        return {int(i): float(np.median(records[records[:, 0] == i, 2]))
                for i in np.unique(records[:, 0]).astype(int)}


def _return_distance(a: np.ndarray, b: np.ndarray) -> float:
    """Symmetric clipped distance; all accepted returns enter the score."""
    if not len(a) and not len(b):
        return 0.0
    if not len(a) or not len(b):
        return 1.0
    # Frequency uses 0.2 MHz and group range uses 25 km. A clipped cost
    # keeps an occasional missed homing frequency from dominating a ridge.
    a_xy = np.column_stack((a[:, 0] / 2.0, a[:, 2] / 25.0))
    b_xy = np.column_stack((b[:, 0] / 2.0, b[:, 2] / 25.0))
    d_ab = cKDTree(b_xy).query(a_xy)[0]
    d_ba = cKDTree(a_xy).query(b_xy)[0]
    return float((np.mean(np.minimum(d_ab / 3.0, 1.0))
                  + np.mean(np.minimum(d_ba / 3.0, 1.0))) / 2.0)


def score(observed: Ionogram, predicted: Ionogram) -> dict[str, float]:
    if not np.array_equal(observed.frequencies, predicted.frequencies):
        raise ValueError("Frequency axes differ")
    if observed.settings != predicted.settings:
        raise ValueError("Ionogram generator settings differ")
    distances = []
    noses = []
    ridges = []
    for mode in (1, -1):
        obs = observed.records[observed.records[:, 1] == mode]
        pred = predicted.records[predicted.records[:, 1] == mode]
        distances.append(_return_distance(obs, pred))
        on, pn = observed.nose(mode), predicted.nose(mode)
        if on is None or pn is None:
            noses.append(1.0)
        elif on >= observed.frequencies[-1] or pn >= observed.frequencies[-1]:
            # A return at the sweep edge gives only a lower bound on the nose.
            noses.append(min(abs(min(on, 10.0) - min(pn, 10.0)) / .8, 1.0)
                         if (on < 10.0 or pn < 10.0) else 0.0)
        else:
            noses.append(min(abs(on - pn) / .8, 1.0))
        oridge, pridge = observed.ridge(mode), predicted.ridge(mode)
        common = sorted(set(oridge) & set(pridge))
        if common:
            residual = np.array([oridge[i] - pridge[i] for i in common])
            ridges.append(float(np.mean(np.minimum(np.abs(residual) / 60.0, 1.0))))
        else:
            ridges.append(0.0 if not oridge and not pridge else 1.0)
    return {
        "return_distance": float(np.mean(distances)),
        "nose": float(np.mean(noses)),
        "ridge": float(np.mean(ridges)),
        "total": float(.55 * np.mean(distances) + .25 * np.mean(noses)
                       + .20 * np.mean(ridges)),
    }


def estimate_density_from_nose(observed: Ionogram, library: list[Ionogram]) -> float:
    """Calibrate the critical-frequency proxy against this exact D generator."""
    estimates = []
    for mode in (1, -1):
        nose = observed.nose(mode)
        pairs = sorted((candidate.nose(mode), candidate.density_scale)
                       for candidate in library
                       if candidate.nose(mode) is not None and candidate.density_scale is not None)
        if nose is None or nose >= 10.0 or len(pairs) < 2:
            continue
        # Average density for each quantized 100 kHz nose, then interpolate.
        unique = sorted(set(item[0] for item in pairs))
        scales = [np.mean([s for n, s in pairs if n == value]) for value in unique]
        estimates.append(float(np.interp(nose, unique, scales)))
    if estimates:
        return float(np.mean(estimates))
    return float(np.median([x.density_scale for x in library if x.density_scale is not None]))


def estimate_height_shift(observed: Ionogram, reference: Ionogram) -> float:
    """Geometric first guess: a 1 km reflector shift changes two-way range ~2 km."""
    differences = []
    for mode in (1, -1):
        obs, ref = observed.ridge(mode), reference.ridge(mode)
        for index in sorted(set(obs) & set(ref)):
            if index <= 40:  # 2–6 MHz: well below the frequency nose
                differences.append(obs[index] - ref[index])
    if len(differences) < 8:
        return 0.0
    return float(np.clip(-np.median(differences) / 2.0, *BOUNDS[1]))


def _record(path: Path, observed: Ionogram, expected: dict | None = None) -> dict:
    ionogram = Ionogram.read(path)
    if ionogram.density_scale is None:
        raise ValueError(f"Missing density_scale in {path}")
    if expected is not None and (abs(ionogram.density_scale - expected["density_scale"]) > 1e-8
                                 or abs(ionogram.hmf2_shift_km - expected["hmf2_shift_km"]) > 1e-8):
        raise ValueError(f"Candidate parameters do not match {path}")
    return {"path": str(path.resolve()), "density_scale": ionogram.density_scale,
            "hmf2_shift_km": ionogram.hmf2_shift_km, "score": score(observed, ionogram)}


def _save(workdir: Path, state: dict) -> None:
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "state.json").write_text(json.dumps(state, indent=2) + "\n")


def _propose(workdir: Path, state: dict) -> None:
    round_number = state["next_round"]
    if round_number > REFINEMENTS:
        return
    rng = np.random.default_rng(20260925 + round_number)
    ranked = sorted(state["evaluations"], key=lambda item: item["score"]["total"])
    best = ranked[0]
    best_xy = np.array([best["density_scale"], best["hmf2_shift_km"]])
    nose_density = state["nose_density_guess"]
    height_guess = state["height_guess"]
    candidates = []
    if round_number == 1:
        center = np.array([nose_density, height_guess])
        span = np.array([.12, 24.0])
        # Stratified coverage about the measured nose and approximate height.
        for i in range(32):
            u = ((i % 8) + rng.uniform()) / 8.0
            v = ((i // 8) + rng.uniform()) / 4.0
            candidates.append(center + span * np.array([2*u-1, 2*v-1]))
        for i in range(8):
            candidates.append(best_xy + rng.normal(size=2) * np.array([.05, 10.0]))
    else:
        # Shrink around the best fits, but retain eight wider candidates to
        # recover from a misleading nose or a local minimum.
        scale = .55 if round_number == 2 else .28
        elite = ranked[:min(8, len(ranked))]
        elite_xy = np.array([[e["density_scale"], e["hmf2_shift_km"]] for e in elite])
        spread = np.maximum(np.std(elite_xy, axis=0), np.array([.025, 4.0]))
        for i in range(32):
            parent = elite_xy[i % len(elite_xy)] if i % 4 == 0 else best_xy
            candidates.append(parent + rng.normal(size=2) * spread * scale)
        for i in range(8):
            candidates.append(best_xy + rng.normal(size=2) * np.array([.065, 14.0]) * scale)
    previous = np.array([[item["density_scale"], item["hmf2_shift_km"]]
                         for item in state["evaluations"]])
    result = []
    for i, point in enumerate(candidates):
        point = np.clip(point, BOUNDS[:, 0], BOUNDS[:, 1])
        # Never spend a costly ray trace on a previous parameter pair.
        attempts = 0
        while np.any(np.all(np.abs(previous - point) < np.array([1e-5, 1e-3]), axis=1)) or any(
                abs(item["density_scale"] - point[0]) < 1e-5
                and abs(item["hmf2_shift_km"] - point[1]) < 1e-3 for item in result):
            point = np.clip(point + rng.normal(size=2) * np.array([.003, .5]),
                            BOUNDS[:, 0], BOUNDS[:, 1])
            attempts += 1
            if attempts > 100:
                raise RuntimeError("Could not propose 40 unique candidates")
        result.append({"task": i + 1, "density_scale": float(point[0]),
                       "hmf2_shift_km": float(point[1])})
    (workdir / f"candidates_round_{round_number}.json").write_text(
        json.dumps(result, indent=2) + "\n")
    print(json.dumps({"round": round_number, "candidate_count": len(result),
                      "best": best, "nose_density_guess": nose_density,
                      "height_guess": height_guess}, indent=2))


def initialize(observed_path: Path, library_dir: Path, workdir: Path) -> None:
    if (workdir / "state.json").exists():
        raise FileExistsError(f"Existing fit in {workdir}")
    observed = Ionogram.read(observed_path)
    paths = sorted(library_dir.glob("*.npz"))
    if len(paths) != POPULATION:
        raise ValueError(f"Expected 40 library ionograms, found {len(paths)}")
    library = [Ionogram.read(path) for path in paths]
    if any(item.hmf2_shift_km != 0.0 for item in library):
        raise ValueError("Initial density library must have zero height shift")
    evaluations = [_record(path, observed) for path in paths]
    density_guess = estimate_density_from_nose(observed, library)
    reference = min(library, key=lambda item: abs(item.density_scale - density_guess))
    state = {"observed": str(observed_path.resolve()), "next_round": 1,
             "nose_density_guess": density_guess,
             "height_guess": estimate_height_shift(observed, reference),
             "evaluations": evaluations}
    _save(workdir, state)
    _propose(workdir, state)


def advance(workdir: Path, results_dir: Path) -> None:
    state = json.loads((workdir / "state.json").read_text())
    round_number = state["next_round"]
    if round_number > REFINEMENTS:
        raise ValueError("All refinement rounds are complete")
    candidates = json.loads((workdir / f"candidates_round_{round_number}.json").read_text())
    observed = Ionogram.read(Path(state["observed"]))
    incoming = []
    for candidate in candidates:
        path = results_dir / f"ionogram_{candidate['task']:02d}.npz"
        incoming.append(_record(path, observed, candidate))
    state["evaluations"].extend(incoming)
    state["next_round"] += 1
    _save(workdir, state)
    _propose(workdir, state)
    if round_number == REFINEMENTS:
        best = min(state["evaluations"], key=lambda item: item["score"]["total"])
        print(json.dumps({"complete": True, "evaluations": len(state["evaluations"]),
                          "best": best}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    init = sub.add_parser("init")
    init.add_argument("observed", type=Path)
    init.add_argument("library_dir", type=Path)
    init.add_argument("workdir", type=Path)
    step = sub.add_parser("advance")
    step.add_argument("workdir", type=Path)
    step.add_argument("results_dir", type=Path)
    args = parser.parse_args()
    if args.action == "init":
        initialize(args.observed, args.library_dir, args.workdir)
    else:
        advance(args.workdir, args.results_dir)


if __name__ == "__main__":
    main()
