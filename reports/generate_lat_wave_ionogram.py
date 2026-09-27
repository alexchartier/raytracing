"""Trace one vertical ionogram from the prebuilt latitude-wave forward grid."""

from __future__ import annotations

import argparse
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
from python_raytrace.multisat_topside_inverse_demo import (  # noqa: E402
    TopsideInverseConfig, _equal_area_vertical_fan, home_frequency_sweep_adaptive,
)

FREQUENCIES = np.arange(2.0, 10.0001, 0.1)


def make_fan() -> tuple[TopsideInverseConfig, np.ndarray, np.ndarray]:
    config = replace(
        TopsideInverseConfig(),
        vertical_elevation_count=24,
        vertical_azimuth_step_deg=30.0,
        vertical_fan_layout="equal_area_guarded",
        vertical_outer_ray_fraction=0.5,
        vertical_guard_seed_limit=4,
        seed_max_candidates_per_frequency=32,
        homed_max_returns_per_frequency=64,
        d_region_model="none",
    )
    elevations, bearings = _equal_area_vertical_fan(config, guard_nadir_rows=3)
    return config, elevations, bearings


def generate(grid_path: Path, latitude: float, longitude: float, altitude: float,
             profile_index: int, output: Path) -> None:
    start = time.perf_counter()
    grid = load_ionosphere_grid_netcdf(grid_path)
    if not (grid.latitudes_deg[0] < latitude < grid.latitudes_deg[-1]
            and grid.longitudes_deg[0] < longitude < grid.longitudes_deg[-1]
            and grid.altitudes_km[0] < altitude < grid.altitudes_km[-1]):
        raise ValueError("Sounder must lie inside the prebuilt forward grid")
    config, elevations, bearings = make_fan()
    rows = np.unique(elevations)
    retained = np.zeros(rows.size, dtype=bool)
    retained[:2] = True
    retained[2::2] = True
    retained[-1] = True
    anchor_directions = int(np.count_nonzero(np.isin(elevations, rows[retained])))
    point = GeoPoint(latitude, longitude, altitude)
    counts = np.zeros((FREQUENCIES.size, 2), dtype=int)
    records = []
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
            records.extend((float(frequency_index), float(mode), ray.group_range_km,
                            ray.miss_m, ray.absorption_db) for ray in rays)
        print(f"profile {profile_index:02d}: mode {mode} complete", flush=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        records=np.asarray(records, dtype=float).reshape(-1, 5),
        count_array=counts,
        frequencies_mhz=FREQUENCIES,
        method=np.array("adaptive"),
        vertical_fan_layout=np.array("equal_area_guarded"),
        vertical_outer_ray_fraction=np.array(0.5),
        vertical_guard_seed_limit=np.array(4),
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
        density_source=np.array("IRI-2016 with imposed wave"),
        runtime_seconds=np.array(time.perf_counter() - start),
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
    args = parser.parse_args()
    generate(args.grid, args.latitude_deg, args.longitude_deg, args.altitude_km,
             args.profile_index, args.output)


if __name__ == "__main__":
    main()
