"""Use a full option-D fan to fill empty, physically plausible pass bins.

The adaptive sweep remains the source of all its accepted returns. A dense
full-fan search is performed only at frequencies where that mode has no
accepted return and the frequency is at most 0.6 MHz above the local peak
plasma frequency. Only reflected returns with group range >=150 km are added.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from generate_lat_wave_ionogram import make_fan, SPACECRAFT_SPEED_MPS, SPACECRAFT_TRACK_BEARING_DEG
from python_raytrace.geometry import GeoPoint
from python_raytrace.grid import load_ionosphere_grid_netcdf
from python_raytrace.multisat_topside_inverse_demo import _home_frequency_returns
from python_raytrace.spacecraft_doppler import spacecraft_doppler
from python_raytrace.tracer import PointToPointRayTracer
from local_ray_lock import local_ray_lock


def recover(source: Path, grid_path: Path, output: Path) -> None:
    start = time.perf_counter()
    with np.load(source, allow_pickle=False) as data:
        payload = {key: data[key] for key in data.files}
    if str(payload["method"]) != "adaptive":
        raise ValueError("Expected the initial adaptive pass ionogram")
    grid = load_ionosphere_grid_netcdf(grid_path)
    latitude = float(payload["tx_lat_deg"])
    longitude = float(payload["tx_lon_deg"])
    altitude = float(payload["tx_alt_km"])
    profile = RegularGridInterpolator(
        (grid.latitudes_deg, grid.longitudes_deg), grid.iono_en_grid,
    )([[latitude, longitude]])[0]
    local_peak_mhz = float(0.00898 * np.sqrt(np.max(profile)))
    frequency_limit = min(10.0, local_peak_mhz + 0.6)
    point = GeoPoint(latitude, longitude, altitude)
    fan_option = (str(payload["vertical_fan_option"])
                  if "vertical_fan_option" in payload else
                  "A" if str(payload["vertical_fan_layout"]) == "az_el" else "D")
    config, elevations, bearings = make_fan(fan_option)
    frequencies = np.asarray(payload["frequencies_mhz"], dtype=float)
    counts = np.asarray(payload["count_array"], dtype=int).copy()
    records = np.asarray(payload["records"], dtype=float).tolist()
    has_spacecraft_doppler = "spacecraft_doppler_hz" in payload
    if has_spacecraft_doppler:
        dopplers = np.asarray(payload["spacecraft_doppler_hz"], dtype=float).tolist()
        launch_angles = np.asarray(payload["launch_angles_deg"], dtype=float).tolist()
        arrival_angles = np.asarray(payload["arrival_angles_deg"], dtype=float).tolist()
        receiver_misses = np.asarray(
            payload["spacecraft_doppler_receiver_miss_m"], dtype=float).tolist()
    tracer = PointToPointRayTracer()
    checks = 0
    recovered_bins = 0
    recovered_returns = 0
    for mode_index, mode in enumerate((1, -1)):
        for index, frequency in enumerate(frequencies):
            if frequency > frequency_limit + 1e-9 or counts[index, mode_index]:
                continue
            checks += 1
            rays = _home_frequency_returns(
                tracer, tx=point, rx=point, grid=grid,
                fan_elevations_deg=elevations, fan_bearings_deg=bearings,
                frequency_mhz=float(frequency), ox_mode=mode, config=config,
                optimizer_method="Powell",
            )
            reflected = [ray for ray in rays if ray.group_range_km >= 150.0]
            if reflected:
                recovered_bins += 1
                recovered_returns += len(reflected)
                counts[index, mode_index] = len(reflected)
                for ray in reflected:
                    records.append([float(index), float(mode), ray.group_range_km,
                                    ray.miss_m, ray.absorption_db])
                    if has_spacecraft_doppler:
                        observable = spacecraft_doppler(
                            ray.ray, point, point, float(frequency),
                            speed_mps=SPACECRAFT_SPEED_MPS,
                            track_bearing_deg=SPACECRAFT_TRACK_BEARING_DEG,
                        )
                        dopplers.append(observable.doppler_hz)
                        launch_angles.append([observable.launch_elevation_deg,
                                              observable.launch_bearing_deg])
                        arrival_angles.append([observable.arrival_elevation_deg,
                                               observable.arrival_bearing_deg])
                        receiver_misses.append(observable.receiver_miss_m)
    payload.update(
        records=np.asarray(records, dtype=float).reshape(-1, 5),
        count_array=counts,
        method=np.array("adaptive_with_dense_gap_recovery"),
        dense_recovery_frequency_limit_mhz=np.array(frequency_limit),
        dense_recovery_frequency_checks=np.array(checks),
        dense_recovered_frequency_bins=np.array(recovered_bins),
        dense_recovered_return_count=np.array(recovered_returns),
        dense_recovery_runtime_seconds=np.array(time.perf_counter() - start),
    )
    if has_spacecraft_doppler:
        payload.update(
            spacecraft_doppler_hz=np.asarray(dopplers, dtype=float),
            launch_angles_deg=np.asarray(launch_angles, dtype=float).reshape(-1, 2),
            arrival_angles_deg=np.asarray(arrival_angles, dtype=float).reshape(-1, 2),
            spacecraft_doppler_receiver_miss_m=np.asarray(receiver_misses, dtype=float),
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **payload)
    print(f"{output}: checked {checks} gaps, recovered {recovered_returns} returns "
          f"in {recovered_bins} bins", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with local_ray_lock():
        recover(args.source, args.grid, args.output)


if __name__ == "__main__":
    main()
