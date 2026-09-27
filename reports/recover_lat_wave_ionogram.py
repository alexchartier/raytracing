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

from generate_lat_wave_ionogram import make_fan
from python_raytrace.geometry import GeoPoint
from python_raytrace.grid import load_ionosphere_grid_netcdf
from python_raytrace.multisat_topside_inverse_demo import _home_frequency_returns
from python_raytrace.tracer import PointToPointRayTracer


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
    config, elevations, bearings = make_fan()
    frequencies = np.asarray(payload["frequencies_mhz"], dtype=float)
    counts = np.asarray(payload["count_array"], dtype=int).copy()
    records = np.asarray(payload["records"], dtype=float).tolist()
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
                records.extend([float(index), float(mode), ray.group_range_km,
                                ray.miss_m, ray.absorption_db] for ray in reflected)
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
    recover(args.source, args.grid, args.output)


if __name__ == "__main__":
    main()
