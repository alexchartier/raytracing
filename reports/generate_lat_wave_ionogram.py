"""Trace one vertical ionogram from the prebuilt latitude-wave forward grid."""

from __future__ import annotations

import argparse
import resource
import sys
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from python_raytrace.geometry import GeoPoint  # noqa: E402
from python_raytrace.grid import load_ionosphere_grid_netcdf  # noqa: E402
from python_raytrace.spacecraft_doppler import spacecraft_doppler  # noqa: E402
from python_raytrace.multisat_topside_inverse_demo import (  # noqa: E402
    TopsideInverseConfig, _equal_area_vertical_fan, _fan_mesh,
    _vertical_search_arrays, home_frequency_sweep_adaptive,
)
from local_ray_lock import local_ray_lock, require_remote_for_large_grid  # noqa: E402

FREQUENCIES = np.arange(2.0, 10.0001, 0.1)
SPACECRAFT_SPEED_MPS = 8000.0
SPACECRAFT_TRACK_BEARING_DEG = 0.0


def make_fan(option: str = "D") -> tuple[TopsideInverseConfig, np.ndarray, np.ndarray]:
    if option not in ("A", "D"):
        raise ValueError("fan option must be A or D")
    config = replace(
        TopsideInverseConfig(),
        vertical_elevation_count=24,
        vertical_azimuth_step_deg=30.0,
        vertical_fan_layout="az_el" if option == "A" else "equal_area_guarded",
        vertical_outer_ray_fraction=1.0 if option == "A" else 0.5,
        vertical_guard_seed_limit=None if option == "A" else 4,
        seed_max_candidates_per_frequency=32,
        homed_max_returns_per_frequency=64,
        d_region_model="none",
    )
    if option == "A":
        elevations, bearings = _fan_mesh(*_vertical_search_arrays(config))
    else:
        elevations, bearings = _equal_area_vertical_fan(config, guard_nadir_rows=3)
    return config, elevations, bearings


def generate(grid_path: Path, latitude: float, longitude: float, altitude: float,
             profile_index: int, output: Path,
             density_source: str = "IRI-2016 with imposed wave",
             fan_option: str = "D") -> None:
    start = time.perf_counter()
    require_remote_for_large_grid(grid_path)
    grid = load_ionosphere_grid_netcdf(grid_path)
    if not (grid.latitudes_deg[0] < latitude < grid.latitudes_deg[-1]
            and grid.longitudes_deg[0] < longitude < grid.longitudes_deg[-1]
            and grid.altitudes_km[0] < altitude < grid.altitudes_km[-1]):
        raise ValueError("Sounder must lie inside the prebuilt forward grid")
    config, elevations, bearings = make_fan(fan_option)
    rows = np.unique(elevations)
    retained = np.zeros(rows.size, dtype=bool)
    retained[:2] = True
    retained[2::2] = True
    retained[-1] = True
    anchor_directions = int(np.count_nonzero(np.isin(elevations, rows[retained])))
    point = GeoPoint(latitude, longitude, altitude)
    counts = np.zeros((FREQUENCIES.size, 2), dtype=int)
    records = []
    dopplers = []
    launch_angles = []
    arrival_angles = []
    receiver_misses = []
    for mode_index, mode in enumerate((1, -1)):
        returns = home_frequency_sweep_adaptive(
            tx=point, rx=point, grid=grid,
            fan_elevations_deg=elevations, fan_bearings_deg=bearings,
            frequencies_mhz=FREQUENCIES, ox_mode=mode, config=config,
            range_min_km=150.0, anchor_stride=5, anchor_block_size=10,
            anchor_elevation_stride=2, anchor_optimizer="Powell",
        )
        for frequency_index, rays in enumerate(returns):
            counts[frequency_index, mode_index] = len(rays)
            for ray in rays:
                observable = spacecraft_doppler(
                    ray.ray, point, point, float(FREQUENCIES[frequency_index]),
                    speed_mps=SPACECRAFT_SPEED_MPS,
                    track_bearing_deg=SPACECRAFT_TRACK_BEARING_DEG,
                )
                records.append((float(frequency_index), float(mode), ray.group_range_km,
                                ray.miss_m, ray.absorption_db))
                dopplers.append(observable.doppler_hz)
                launch_angles.append((observable.launch_elevation_deg,
                                      observable.launch_bearing_deg))
                arrival_angles.append((observable.arrival_elevation_deg,
                                       observable.arrival_bearing_deg))
                receiver_misses.append(observable.receiver_miss_m)
        print(f"profile {profile_index:02d}: mode {mode} complete", flush=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        records=np.asarray(records, dtype=float).reshape(-1, 5),
        spacecraft_doppler_hz=np.asarray(dopplers, dtype=float),
        launch_angles_deg=np.asarray(launch_angles, dtype=float).reshape(-1, 2),
        arrival_angles_deg=np.asarray(arrival_angles, dtype=float).reshape(-1, 2),
        spacecraft_doppler_receiver_miss_m=np.asarray(receiver_misses, dtype=float),
        spacecraft_speed_mps=np.array(SPACECRAFT_SPEED_MPS),
        spacecraft_track_bearing_deg=np.array(SPACECRAFT_TRACK_BEARING_DEG),
        doppler_model=np.array("f/c times velocity dot launch-minus-arrival phase momentum"),
        count_array=counts,
        frequencies_mhz=FREQUENCIES,
        method=np.array("adaptive"),
        vertical_fan_option=np.array(fan_option),
        vertical_fan_layout=np.array(config.vertical_fan_layout),
        vertical_outer_ray_fraction=np.array(config.vertical_outer_ray_fraction),
        vertical_guard_seed_limit=np.array(-1 if config.vertical_guard_seed_limit is None
                                           else config.vertical_guard_seed_limit),
        anchor_stride=np.array(5),
        anchor_block_size=np.array(10),
        anchor_elevation_stride=np.array(2),
        anchor_optimizer=np.array("Powell"),
        fan_launch_directions=np.array(len(elevations)),
        anchor_fan_launch_directions=np.array(anchor_directions),
        homing_tolerance_m=np.array(config.homing_tolerance_m),
        tx_lat_deg=np.array(latitude),
        tx_lon_deg=np.array(longitude),
        tx_alt_km=np.array(altitude),
        profile_index=np.array(profile_index),
        density_source=np.array(density_source),
        runtime_seconds=np.array(time.perf_counter() - start),
        peak_rss_mb=np.array(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576.0),
    )
    print(f"{output}: {len(records)} accepted returns", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--latitude-deg", type=float, required=True)
    parser.add_argument("--longitude-deg", type=float, required=True)
    parser.add_argument("--altitude-km", type=float, required=True)
    parser.add_argument("--profile-index", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--density-source", default="IRI-2016 with imposed wave")
    parser.add_argument("--fan-option", choices=("A", "D"), default="D")
    args = parser.parse_args()
    with local_ray_lock():
        generate(args.grid, args.latitude_deg, args.longitude_deg, args.altitude_km,
                 args.profile_index, args.output, args.density_source, args.fan_option)


if __name__ == "__main__":
    main()
