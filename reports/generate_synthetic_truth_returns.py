"""Trace the dense-fan synthetic truth sweep used by the standalone ionogram.

Run batches with ``python3 reports/generate_synthetic_truth_returns.py --start 0 --stop 10``
and then use ``--merge`` after all 0–10, 10–20, ..., 80–81 batches finish.
Requires the installed PyIRI/PHaRLAP runtime.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reports.build_topside_validation_report import write_orbit_fixture  # noqa: E402
from python_raytrace.multisat_topside_inverse_demo import (  # noqa: E402
    IonosphereFitParams,
    TopsideInverseConfig,
    _apply_fit_params_to_grid,
    _home_frequency_returns,
    _subset_problem_cases,
    build_inverse_problem,
)
from python_raytrace.tracer import PointToPointRayTracer  # noqa: E402


OUTPUT = ROOT / "reports" / "data" / "synthetic_truth_returns_dense_2-10MHz_100kHz.npz"
FREQUENCIES = np.arange(2.0, 10.0001, 0.1)


def chunk_path(start: int, stop: int) -> Path:
    return OUTPUT.parent / f"synthetic_truth_dense_chunk_{start:02d}_{stop:02d}.npz"


def merge_chunks() -> None:
    chunk_ranges = [(start, min(start + 10, len(FREQUENCIES))) for start in range(0, len(FREQUENCIES), 10)]
    record_chunks = []
    counts = np.zeros((len(FREQUENCIES), 2), dtype=int)
    for start, stop in chunk_ranges:
        with np.load(chunk_path(start, stop), allow_pickle=False) as chunk:
            record_chunks.append(np.asarray(chunk["records"], dtype=float))
            counts += np.asarray(chunk["count_array"], dtype=int)
    records = np.concatenate(record_chunks)
    np.savez_compressed(
        OUTPUT,
        records=records,
        count_array=counts,
        frequencies_mhz=FREQUENCIES,
        fan_launch_directions=np.array(288),
        homing_tolerance_m=np.array(1000.0),
    )
    print(f"accepted: {len(records)}; at or above 150 km: "
          f"{np.count_nonzero(records[:, 2] >= 150.0)}", flush=True)
    print(OUTPUT, flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--merge", action="store_true")
    args = parser.parse_args()
    if args.merge:
        merge_chunks()
        return
    if args.start is None or args.stop is None or not 0 <= args.start < args.stop <= len(FREQUENCIES):
        parser.error("choose a batch with --start and --stop between 0 and 81")
    with tempfile.TemporaryDirectory(prefix="topside_truth_sweep_") as directory:
        temporary = Path(directory)
        orbit_path = temporary / "generated_orbit.nc"
        write_orbit_fixture(orbit_path)
        config = replace(
            TopsideInverseConfig(),
            ampere_file=orbit_path,
            planes=(1,),
            frequencies_mhz=(4.0, 5.0, 6.0),
            vertical_elevation_count=24,
            vertical_azimuth_step_deg=30.0,
            oblique_elevation_count=9,
            oblique_bearing_count=7,
            seed_max_candidates_per_frequency=32,
            homed_max_returns_per_frequency=64,
            grid_lat_step_deg=2.0,
            grid_lon_step_deg=4.0,
            grid_alt_step_km=20.0,
            d_region_model="none",
            grid_cache_path=temporary / "background.nc",
        )
        problem = _subset_problem_cases(build_inverse_problem(config), ["plane01_sv121_vertical"])
        case = problem.cases[0]
        grid = _apply_fit_params_to_grid(
            problem,
            problem.background_grids[0],
            IonosphereFitParams(
                density_scale=1.12,
                hmf2_shift_km=0.0,
                wave_amplitude_fraction=0.0,
                wave_phase_rad=0.0,
                wave_bearing_deg=0.0,
            ),
        )
        counts = np.zeros((FREQUENCIES.size, 2), dtype=int)
        records: list[tuple[float, float, float, float, float]] = []
        for frequency_index in range(args.start, args.stop):
            frequency_mhz = FREQUENCIES[frequency_index]
            tracer = PointToPointRayTracer()
            for mode_index, mode in enumerate((1, -1)):
                returns = _home_frequency_returns(
                    tracer,
                    tx=case.tx_points[1],
                    rx=case.rx_points[1],
                    grid=grid,
                    fan_elevations_deg=case.fan_elevations_deg,
                    fan_bearings_deg=case.fan_bearings_deg,
                    frequency_mhz=float(frequency_mhz),
                    ox_mode=mode,
                    config=config,
                )
                counts[frequency_index, mode_index] = len(returns)
                records.extend(
                    (float(frequency_index), float(mode), ray.group_range_km,
                     ray.miss_m, ray.absorption_db)
                    for ray in returns
                )
            print(f"{frequency_index + 1}/{FREQUENCIES.size} frequencies", flush=True)

        records_array = np.asarray(records, dtype=float).reshape(-1, 5)
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            chunk_path(args.start, args.stop),
            records=records_array,
            count_array=counts,
        )
        print(chunk_path(args.start, args.stop), flush=True)


if __name__ == "__main__":
    main()
