"""Trace synthetic truth returns used by the standalone ionogram.

Run adaptive batches with ``python3 reports/generate_synthetic_truth_returns.py
--start 0 --stop 10`` and then use ``--merge`` after all 0–10, 10–20,
..., 80–81 batches finish. ``--method dense`` reproduces the original sweep.
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
    home_frequency_sweep_adaptive,
)
from python_raytrace.tracer import PointToPointRayTracer  # noqa: E402


DATA_DIR = ROOT / "reports" / "data"
FREQUENCIES = np.arange(2.0, 10.0001, 0.1)


def output_path(method: str) -> Path:
    return DATA_DIR / f"synthetic_truth_returns_{method}_2-10MHz_100kHz.npz"


def chunk_path(method: str, start: int, stop: int) -> Path:
    return DATA_DIR / f"synthetic_truth_{method}_chunk_{start:02d}_{stop:02d}.npz"


def merge_chunks(method: str, anchor_stride: int) -> None:
    chunk_ranges = [(start, min(start + 10, len(FREQUENCIES))) for start in range(0, len(FREQUENCIES), 10)]
    record_chunks = []
    counts = np.zeros((len(FREQUENCIES), 2), dtype=int)
    for start, stop in chunk_ranges:
        with np.load(chunk_path(method, start, stop), allow_pickle=False) as chunk:
            if str(chunk["method"]) != method or int(chunk["anchor_stride"]) != anchor_stride:
                raise ValueError(f"Inconsistent sweep settings in {chunk_path(method, start, stop)}")
            record_chunks.append(np.asarray(chunk["records"], dtype=float))
            counts += np.asarray(chunk["count_array"], dtype=int)
    records = np.concatenate(record_chunks)
    output = output_path(method)
    np.savez_compressed(
        output,
        records=records,
        count_array=counts,
        frequencies_mhz=FREQUENCIES,
        method=np.array(method),
        anchor_stride=np.array(anchor_stride),
        fan_launch_directions=np.array(288),
        homing_tolerance_m=np.array(1000.0),
    )
    print(f"accepted: {len(records)}; at or above 150 km: "
          f"{np.count_nonzero(records[:, 2] >= 150.0)}", flush=True)
    print(output, flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--method", choices=("adaptive", "dense"), default="adaptive")
    parser.add_argument("--anchor-stride", type=int, default=5)
    args = parser.parse_args()
    if args.anchor_stride < 1:
        parser.error("--anchor-stride must be positive")
    if args.merge:
        merge_chunks(args.method, args.anchor_stride)
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
        for mode_index, mode in enumerate((1, -1)):
            if args.method == "adaptive":
                return_sets = home_frequency_sweep_adaptive(
                    tx=case.tx_points[1], rx=case.rx_points[1], grid=grid,
                    fan_elevations_deg=case.fan_elevations_deg,
                    fan_bearings_deg=case.fan_bearings_deg,
                    frequencies_mhz=FREQUENCIES[args.start:args.stop],
                    ox_mode=mode, config=config,
                    range_min_km=150.0, anchor_stride=args.anchor_stride,
                )
            else:
                return_sets = tuple(
                    _home_frequency_returns(
                        PointToPointRayTracer(),
                        tx=case.tx_points[1], rx=case.rx_points[1], grid=grid,
                        fan_elevations_deg=case.fan_elevations_deg,
                        fan_bearings_deg=case.fan_bearings_deg,
                        frequency_mhz=float(FREQUENCIES[frequency_index]),
                        ox_mode=mode, config=config,
                    )
                    for frequency_index in range(args.start, args.stop)
                )
            for frequency_index, returns in zip(range(args.start, args.stop), return_sets):
                counts[frequency_index, mode_index] = len(returns)
                records.extend(
                    (float(frequency_index), float(mode), ray.group_range_km,
                     ray.miss_m, ray.absorption_db)
                    for ray in returns
                )
                print(f"{frequency_index + 1}/{FREQUENCIES.size} frequencies, mode {mode}", flush=True)

        records_array = np.asarray(records, dtype=float).reshape(-1, 5)
        DATA_DIR.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            chunk_path(args.method, args.start, args.stop),
            records=records_array,
            count_array=counts,
            method=np.array(args.method),
            anchor_stride=np.array(args.anchor_stride),
        )
        print(chunk_path(args.method, args.start, args.stop), flush=True)


if __name__ == "__main__":
    main()
