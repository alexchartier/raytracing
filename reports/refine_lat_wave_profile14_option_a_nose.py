"""Recover profile-14 option-A O returns by fine frequency continuation.

The 288-direction A fan is first homed at 4.6 and 5.1 MHz. Accepted paths are
then continued in 20 kHz increments. Only 100 kHz grid frequencies and returns
with group range at least 150 km enter the saved ionogram.
"""

from __future__ import annotations

import json
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

DATA = ROOT / "reports/data"
SOURCE = DATA / "lat_wave_profile14_option_a_truth_recovered.npz"
OUTPUT = DATA / "lat_wave_profile14_option_a_truth_nose_refined.npz"
GRID = DATA / "lat_wave_pass_forward_grid.nc"
MANIFEST = DATA / "lat_wave_pass_manifest.json"


def _direction(elevation_deg: float, bearing_deg: float) -> np.ndarray:
    elevation = np.deg2rad(elevation_deg)
    bearing = np.deg2rad(bearing_deg)
    return np.array([np.cos(elevation) * np.sin(bearing),
                     np.cos(elevation) * np.cos(bearing), np.sin(elevation)])


def main() -> None:
    started = time.perf_counter()
    with np.load(SOURCE, allow_pickle=False) as source:
        payload = {key: source[key] for key in source.files}
    if (int(payload["profile_index"]) != 14
            or str(payload["vertical_fan_option"]) != "A"):
        raise ValueError("Expected recovered option-A truth profile 14")
    profile = json.loads(MANIFEST.read_text())["profiles"][13]
    point = GeoPoint(profile["latitude_deg"], profile["longitude_deg"],
                     profile["altitude_km"])
    grid = load_ionosphere_grid_netcdf(GRID)
    config, elevations, bearings = make_fan("A")
    tracer = PointToPointRayTracer()
    found = []
    branch_status = []
    for start_centimhz, stop_centimhz, step in (
            (460, 450, -2), (460, 510, 2),
            (510, 450, -2), (510, 600, 2)):
        start_mhz = start_centimhz / 100.0
        initial = _home_frequency_returns(
            tracer, tx=point, rx=point, grid=grid,
            fan_elevations_deg=elevations, fan_bearings_deg=bearings,
            frequency_mhz=start_mhz, ox_mode=1, config=config,
            optimizer_method="Powell",
        )
        paths = tuple(ray for ray in initial if ray.group_range_km >= 150.0)
        if not paths:
            raise RuntimeError(f"No physical seed return at {start_mhz:.1f} MHz")
        for centimhz in range(start_centimhz + step, stop_centimhz + step, step):
            frequency_mhz = centimhz / 100.0
            paths = _continue_homed_returns(
                tracer, paths, tx=point, rx=point, grid=grid,
                frequency_mhz=frequency_mhz, ox_mode=1, config=config,
                range_min_km=150.0,
            )
            if not paths:
                branch_status.append({"seed_mhz": start_mhz,
                                      "direction": "up" if step > 0 else "down",
                                      "last_attempted_mhz": frequency_mhz,
                                      "continued": False})
                break
            if centimhz % 10 == 0:
                frequency_index = (centimhz - 200) // 10
                found.extend((frequency_index, ray) for ray in paths)
        else:
            branch_status.append({"seed_mhz": start_mhz,
                                  "direction": "up" if step > 0 else "down",
                                  "last_attempted_mhz": stop_centimhz / 100.0,
                                  "continued": True})

    records = np.asarray(payload["records"], dtype=float).tolist()
    dopplers = np.asarray(payload["spacecraft_doppler_hz"], dtype=float).tolist()
    launch_angles = np.asarray(payload["launch_angles_deg"], dtype=float).tolist()
    arrival_angles = np.asarray(payload["arrival_angles_deg"], dtype=float).tolist()
    receiver_misses = np.asarray(payload["spacecraft_doppler_receiver_miss_m"],
                                 dtype=float).tolist()
    counts = np.asarray(payload["count_array"], dtype=int).copy()
    added = []
    for frequency_index, ray in sorted(found, key=lambda item: (item[0], item[1].group_range_km)):
        observable = spacecraft_doppler(ray.ray, point, point,
                                        2.0 + 0.1 * frequency_index)
        direction = _direction(observable.launch_elevation_deg,
                               observable.launch_bearing_deg)
        duplicate = False
        for record, angles in zip(records, launch_angles):
            if (int(record[0]) == frequency_index and int(record[1]) == 1
                    and abs(record[2] - ray.group_range_km) < 1.0
                    and np.dot(direction, _direction(*angles))
                    >= np.cos(np.deg2rad(0.02))):
                duplicate = True
                break
        if duplicate:
            continue
        records.append([float(frequency_index), 1.0, ray.group_range_km,
                        ray.miss_m, ray.absorption_db])
        dopplers.append(observable.doppler_hz)
        launch_angles.append([observable.launch_elevation_deg,
                              observable.launch_bearing_deg])
        arrival_angles.append([observable.arrival_elevation_deg,
                               observable.arrival_bearing_deg])
        receiver_misses.append(observable.receiver_miss_m)
        counts[frequency_index, 0] += 1
        added.append({"frequency_mhz": 2.0 + 0.1 * frequency_index,
                      "group_range_km": ray.group_range_km,
                      "homing_miss_m": ray.miss_m})
    payload.update(
        records=np.asarray(records, dtype=float).reshape(-1, 5),
        spacecraft_doppler_hz=np.asarray(dopplers, dtype=float),
        launch_angles_deg=np.asarray(launch_angles, dtype=float).reshape(-1, 2),
        arrival_angles_deg=np.asarray(arrival_angles, dtype=float).reshape(-1, 2),
        spacecraft_doppler_receiver_miss_m=np.asarray(receiver_misses, dtype=float),
        count_array=counts,
        method=np.array("adaptive_with_dense_gap_recovery_and_20khz_o_nose_continuation"),
        nose_refinement_added_return_count=np.array(len(added)),
        nose_refinement_runtime_seconds=np.array(time.perf_counter() - started),
        nose_refinement_accepted_returns_json=np.array(json.dumps(added)),
        nose_refinement_branch_status_json=np.array(json.dumps(branch_status)),
    )
    np.savez_compressed(OUTPUT, **payload)
    print(json.dumps({"output": str(OUTPUT), "added": added,
                      "branch_status": branch_status,
                      "runtime_seconds": time.perf_counter() - started}, indent=2))


if __name__ == "__main__":
    main()
