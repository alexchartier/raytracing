"""Recover near-nose O/X returns with 20 kHz frequency continuation.

The input is an already traced 100 kHz ionogram. Full-fan seeds are obtained
at its measured nose and 0.5 and 1.0 MHz below it. Each seed is followed up
and down by up to 0.54 MHz in 20 kHz steps. Only rays actually traced and homed at
100 kHz grid frequencies enter the output. The density grid is used only to
raytrace, never to set a cutoff or fit parameter.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from reports.generate_lat_wave_ionogram import make_fan  # noqa: E402
from python_raytrace.geometry import GeoPoint  # noqa: E402
from python_raytrace.grid import load_ionosphere_grid_netcdf  # noqa: E402
from python_raytrace.multisat_topside_inverse_demo import (  # noqa: E402
    _continue_homed_returns, _home_frequency_returns,
)
from python_raytrace.spacecraft_doppler import spacecraft_doppler  # noqa: E402
from python_raytrace.tracer import PointToPointRayTracer  # noqa: E402
from reports.local_ray_lock import local_ray_lock  # noqa: E402


def direction(elevation_deg: float, bearing_deg: float) -> np.ndarray:
    elevation = np.deg2rad(elevation_deg)
    bearing = np.deg2rad(bearing_deg)
    return np.array([np.cos(elevation) * np.sin(bearing),
                     np.cos(elevation) * np.cos(bearing), np.sin(elevation)])


def continue_ionogram(source_path: Path, grid_path: Path, output_path: Path,
                      above_only: bool = False) -> dict:
    started = time.perf_counter()
    with np.load(source_path, allow_pickle=False) as source:
        payload = {key: source[key] for key in source.files}
    frequencies = np.asarray(payload["frequencies_mhz"], dtype=float)
    if len(frequencies) != 81 or not np.allclose(frequencies, 2.0 + 0.1 * np.arange(81)):
        raise ValueError("Expected the 2–10 MHz, 100 kHz grid")
    for key in ("spacecraft_doppler_hz", "launch_angles_deg",
                "arrival_angles_deg", "spacecraft_doppler_receiver_miss_m"):
        if key not in payload:
            raise ValueError(f"Missing {key} from source ionogram")
    source_count = len(payload["records"])
    if (any(len(payload[key]) != source_count for key in
            ("spacecraft_doppler_hz", "launch_angles_deg",
             "arrival_angles_deg", "spacecraft_doppler_receiver_miss_m"))
            or np.asarray(payload["launch_angles_deg"]).shape != (source_count, 2)
            or np.asarray(payload["arrival_angles_deg"]).shape != (source_count, 2)):
        raise ValueError("Source ray records and Doppler/angle arrays do not align")
    if above_only and "20khz_nose_continuation" not in str(payload["method"]):
        raise ValueError("Above-nose probes require a continued source ionogram")
    option = (str(payload["vertical_fan_option"])
              if "vertical_fan_option" in payload else
              "A" if str(payload["vertical_fan_layout"]) == "az_el" else "D")
    config, elevations, bearings = make_fan(option)
    grid = load_ionosphere_grid_netcdf(grid_path)
    point = GeoPoint(float(payload["tx_lat_deg"]), float(payload["tx_lon_deg"]),
                     float(payload["tx_alt_km"]))
    tracer = PointToPointRayTracer()
    records = np.asarray(payload["records"], dtype=float).tolist()
    dopplers = np.asarray(payload["spacecraft_doppler_hz"], dtype=float).tolist()
    launch_angles = np.asarray(payload["launch_angles_deg"], dtype=float).tolist()
    arrival_angles = np.asarray(payload["arrival_angles_deg"], dtype=float).tolist()
    receiver_misses = np.asarray(payload["spacecraft_doppler_receiver_miss_m"],
                                 dtype=float).tolist()
    counts = np.asarray(payload["count_array"], dtype=int).copy()
    found = []
    seeded = []
    branches = []
    for mode in (1, -1):
        mode_records = np.asarray(records, dtype=float)
        mode_records = mode_records[mode_records[:, 1] == mode]
        if len(mode_records) == 0:
            continue
        nose_index = int(max(mode_records[:, 0]))
        # 10 centi-MHz is one 100 kHz ionogram bin; 2 is a 20 kHz step.
        seed_indices = (sorted(set(max(0, nose_index - offset)
                                   for offset in (10, 5, 0))) if not above_only
                        else sorted(set(min(80, nose_index + offset)
                                        for offset in (2, 4, 6))))
        for seed_index in seed_indices:
            seed_centimhz = 200 + 10 * seed_index
            initial = tuple(ray for ray in _home_frequency_returns(
                tracer, tx=point, rx=point, grid=grid,
                fan_elevations_deg=elevations, fan_bearings_deg=bearings,
                frequency_mhz=seed_centimhz / 100.0,
                ox_mode=mode, config=config, optimizer_method="Powell",
            ) if ray.group_range_km >= 150.0)
            seeded.append({"mode": mode, "frequency_mhz": seed_centimhz / 100.0,
                           "accepted": len(initial)})
            found.extend((seed_index, mode, ray) for ray in initial)
            for step in (-2, 2):
                paths = initial
                attempted = seed_centimhz
                for centimhz in range(seed_centimhz + step,
                                      seed_centimhz + step * 28, step):
                    if centimhz < 200 or centimhz > 1000:
                        break
                    attempted = centimhz
                    paths = _continue_homed_returns(
                        tracer, paths, tx=point, rx=point, grid=grid,
                        frequency_mhz=centimhz / 100.0, ox_mode=mode,
                        config=config, range_min_km=150.0,
                    )
                    if not paths:
                        break
                    if centimhz % 10 == 0:
                        frequency_index = (centimhz - 200) // 10
                        found.extend((frequency_index, mode, ray) for ray in paths)
                branches.append({"mode": mode, "seed_mhz": seed_centimhz / 100.0,
                                 "direction": "up" if step > 0 else "down",
                                 "last_attempted_mhz": attempted / 100.0,
                                 "last_step_homed": bool(paths)})
    added = []
    speed_mps = float(payload["spacecraft_speed_mps"])
    bearing_deg = float(payload["spacecraft_track_bearing_deg"])
    for index, mode, ray in sorted(found, key=lambda row:
                                   (row[0], row[1], row[2].group_range_km)):
        frequency_mhz = float(frequencies[index])
        if not math.isclose(ray.frequency_mhz, frequency_mhz, abs_tol=1e-9):
            raise AssertionError("Continuation ray was not traced at its plotted frequency")
        observable = spacecraft_doppler(ray.ray, point, point, frequency_mhz,
                                        speed_mps=speed_mps,
                                        track_bearing_deg=bearing_deg)
        new_direction = direction(observable.launch_elevation_deg,
                                  observable.launch_bearing_deg)
        if any(int(record[0]) == index and int(record[1]) == mode
               and abs(record[2] - ray.group_range_km) < 1.0
               and float(np.dot(new_direction, direction(*angles)))
               >= math.cos(math.radians(0.02))
               for record, angles in zip(records, launch_angles)):
            continue
        records.append([float(index), float(mode), ray.group_range_km,
                        ray.miss_m, ray.absorption_db])
        dopplers.append(observable.doppler_hz)
        launch_angles.append([observable.launch_elevation_deg,
                              observable.launch_bearing_deg])
        arrival_angles.append([observable.arrival_elevation_deg,
                               observable.arrival_bearing_deg])
        receiver_misses.append(observable.receiver_miss_m)
        counts[index, 0 if mode == 1 else 1] += 1
        added.append({"frequency_mhz": frequency_mhz, "mode": mode,
                      "group_range_km": ray.group_range_km,
                      "homing_miss_m": ray.miss_m})
    elapsed = time.perf_counter() - started
    payload.update(
        records=np.asarray(records, dtype=float).reshape(-1, 5),
        spacecraft_doppler_hz=np.asarray(dopplers, dtype=float),
        launch_angles_deg=np.asarray(launch_angles, dtype=float).reshape(-1, 2),
        arrival_angles_deg=np.asarray(arrival_angles, dtype=float).reshape(-1, 2),
        spacecraft_doppler_receiver_miss_m=np.asarray(receiver_misses, dtype=float),
        count_array=counts,
    )
    if above_only:
        payload.update(
            method=np.array("adaptive_with_dense_gap_recovery_and_20khz_nose_continuation_and_above_nose_probes"),
            above_nose_added_return_count=np.array(len(added)),
            above_nose_runtime_seconds=np.array(elapsed),
            above_nose_accepted_returns_json=np.array(json.dumps(added)),
            above_nose_seeds_json=np.array(json.dumps(seeded)),
            above_nose_branches_json=np.array(json.dumps(branches)),
        )
    else:
        payload.update(
            method=np.array("adaptive_with_dense_gap_recovery_and_20khz_nose_continuation"),
            nose_refinement_added_return_count=np.array(len(added)),
            nose_refinement_runtime_seconds=np.array(elapsed),
            nose_refinement_accepted_returns_json=np.array(json.dumps(added)),
            nose_refinement_seeds_json=np.array(json.dumps(seeded)),
            nose_refinement_branches_json=np.array(json.dumps(branches)),
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **payload)
    return {"profile_index": int(payload["profile_index"]), "source": str(source_path),
            "output": str(output_path), "original_returns": len(payload["records"]) - len(added),
            "added_returns": len(added), "added_by_mode": {
                str(mode): sum(row["mode"] == mode for row in added) for mode in (1, -1)},
            "seeded": seeded, "runtime_seconds": elapsed,
            "above_only": above_only}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--above-only", action="store_true",
                        help="Probe 0.2, 0.4, and 0.6 MHz above the source noses")
    args = parser.parse_args()
    with local_ray_lock():
        print(json.dumps(continue_ionogram(args.source, args.grid, args.output,
                                           args.above_only), indent=2))


if __name__ == "__main__":
    main()
