"""Re-trace interior frequency gaps with 20 kHz ray continuation.

``trace`` runs on Cartman with a private density grid. ``merge`` runs locally
and adds only accepted rays from bins empty in the source ionogram. Every
added return was traced at its exact 100 kHz display frequency.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def interior_gaps(records: np.ndarray, frequencies: np.ndarray, mode: int) -> list[list[int]]:
    present = sorted(set(records[records[:, 1] == mode, 0].astype(int)))
    if len(present) < 2:
        return []
    missing = sorted(set(range(present[0] + 1, present[-1])) - set(present))
    groups: list[list[int]] = []
    for index in missing:
        if not groups or index != groups[-1][-1] + 1:
            groups.append([index])
        else:
            groups[-1].append(index)
    return groups


def trace(source: Path, density_grid: Path, output: Path) -> None:
    from reports import generate_synthetic_truth_returns as generator
    from python_raytrace.multisat_topside_inverse_demo import (
        _continue_homed_returns, _deduplicate_homed_returns,
        _home_frequency_returns,
    )
    from python_raytrace.tracer import PointToPointRayTracer

    with np.load(source, allow_pickle=False) as data:
        frequencies = np.asarray(data["frequencies_mhz"], dtype=float)
        original = np.asarray(data["records"], dtype=float)
    if len(frequencies) != 81 or not np.allclose(frequencies, 2 + .1 * np.arange(81)):
        raise ValueError("Expected a complete 2–10 MHz, 100 kHz ionogram")
    plans = {mode: interior_gaps(original, frequencies, mode) for mode in (1, -1)}
    print("gap plans", {m: [[round(float(frequencies[i]), 2) for i in group]
                                for group in groups] for m, groups in plans.items()}, flush=True)
    if not any(plans.values()):
        np.savez_compressed(output, frequencies_mhz=frequencies,
                            records=np.empty((0, 5), dtype=float),
                            count_array=np.zeros((len(frequencies), 2), dtype=int))
        print("no interior gaps; no rays needed", flush=True)
        return

    def recover(**kwargs):
        mode = kwargs["ox_mode"]
        groups = plans[mode]
        results = [[] for _ in frequencies]
        if not groups:
            return tuple(() for _ in frequencies)
        tracer = PointToPointRayTracer(cache_native_grid=True)

        def full_fan(frequency: float):
            return tuple(ray for ray in _home_frequency_returns(
                tracer, tx=kwargs["tx"], rx=kwargs["rx"], grid=kwargs["grid"],
                fan_elevations_deg=kwargs["fan_elevations_deg"],
                fan_bearings_deg=kwargs["fan_bearings_deg"],
                frequency_mhz=frequency, ox_mode=mode, config=kwargs["config"],
                optimizer_method="Powell") if ray.group_range_km >= 150)

        for group in groups:
            start, stop = group[0], group[-1]
            for seed_index, direction in ((start - 1, 1), (stop + 1, -1)):
                seed_frequency = round(float(frequencies[seed_index]), 2)
                paths = full_fan(seed_frequency)
                print("seed", mode, seed_frequency, len(paths), flush=True)
                for step in range(1, 10 * (len(group) + 1) + 1):
                    frequency = round(seed_frequency + direction * .02 * step, 2)
                    paths = _continue_homed_returns(
                        tracer, paths, tx=kwargs["tx"], rx=kwargs["rx"],
                        grid=kwargs["grid"], frequency_mhz=frequency,
                        ox_mode=mode, config=kwargs["config"], range_min_km=150)
                    if not paths:
                        break
                    match = np.flatnonzero(np.isclose(frequencies[group], frequency,
                                                       rtol=0, atol=1e-8))
                    for local in match:
                        results[group[int(local)]].extend(paths)
                    if (direction > 0 and frequency >= frequencies[stop] - 1e-8 or
                            direction < 0 and frequency <= frequencies[start] + 1e-8):
                        break
            # A dense fan at the exact display bin is the fallback if both
            # neighbouring continuations failed.
            for index in group:
                if not results[index]:
                    results[index].extend(full_fan(float(frequencies[index])))
            print("gap result", mode, [(round(float(frequencies[i]), 2),
                                        len(results[i])) for i in group], flush=True)
        return tuple(_deduplicate_homed_returns(
            rays, kwargs["config"].homed_max_returns_per_frequency)
                     for rays in results)

    generator.home_frequency_sweep_adaptive = recover
    sys.argv = ["generate_synthetic_truth_returns.py", "--start", "0", "--stop", "81",
                "--output", str(output), "--density-grid-npz", str(density_grid),
                "--density-scale", "1.0"]
    generator.main()


def merge(source: Path, recovered: Path, output: Path, strict: bool) -> None:
    with np.load(source, allow_pickle=False) as data:
        payload = {key: data[key] for key in data.files}
    with np.load(recovered, allow_pickle=False) as data:
        supplement = {key: data[key] for key in data.files}
    frequencies = np.asarray(payload["frequencies_mhz"], dtype=float)
    if not np.allclose(frequencies, supplement["frequencies_mhz"], rtol=0, atol=1e-10):
        raise ValueError("Recovered rays have a different frequency axis")
    original = np.asarray(payload["records"], dtype=float)
    added = np.asarray(supplement["records"], dtype=float)
    allowed = {(index, mode) for mode in (1, -1)
               for group in interior_gaps(original, frequencies, mode)
               for index in group}
    gate = float(payload["homing_tolerance_m"])
    if (np.any(~np.isfinite(added)) or np.any(added[:, 2] < 150)
            or np.any(added[:, 3] > gate)
            or any((int(row[0]), int(row[1])) not in allowed for row in added)):
        raise ValueError("Recovered ray fails the original gap or homing gate")
    found = {(int(row[0]), int(row[1])) for row in added}
    unrecovered = sorted(allowed - found)
    if strict and unrecovered:
        raise RuntimeError(f"Interior gaps remain at {unrecovered}")
    payload["records"] = np.asarray(sorted(np.concatenate((original, added)),
                                          key=lambda row: (row[0], row[1], row[2])))
    counts = np.asarray(payload["count_array"], dtype=int).copy()
    for row in added:
        counts[int(row[0]), 0 if int(row[1]) == 1 else 1] += 1
    payload["count_array"] = counts
    payload["method"] = np.array("adaptive_with_dense_gap_recovery")
    provenance = {"source": source.name, "recovered": recovered.name,
                  "accepted_added": len(added), "unrecovered_bins": unrecovered,
                  "frequency_step_mhz": .02, "homing_gate_m": gate}
    payload["gap_recovery_json"] = np.array(json.dumps(provenance))
    np.savez_compressed(output, **payload)
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    tracing = sub.add_parser("trace")
    tracing.add_argument("--source", type=Path, required=True)
    tracing.add_argument("--density-grid", type=Path, required=True)
    tracing.add_argument("--output", type=Path, required=True)
    merging = sub.add_parser("merge")
    merging.add_argument("--source", type=Path, required=True)
    merging.add_argument("--recovered", type=Path, required=True)
    merging.add_argument("--output", type=Path, required=True)
    merging.add_argument("--strict", action="store_true")
    args = parser.parse_args()
    if args.action == "trace":
        trace(args.source, args.density_grid, args.output)
    else:
        merge(args.source, args.recovered, args.output, args.strict)
