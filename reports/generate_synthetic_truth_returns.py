"""Trace synthetic truth returns used by the standalone ionogram.

Run the complete adaptive sweep with ``python3
reports/generate_synthetic_truth_returns.py --start 0 --stop 81``. The default
vertical fan is option D: half-density equal-area cells outside three near-nadir
guard rings, with up to four extra guard seeds. Its output uses a distinct ``_D``
suffix. Use ``--vertical-fan-layout az_el`` for the original rectangular fan. Shorter
batches can be merged with ``--merge``. ``--method dense`` traces the full fan.
Requires the installed PyIRI/PHaRLAP runtime.
"""

from __future__ import annotations

import argparse
import sys
import tempfile
import time
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


def output_path(method: str, fan_layout: str = "az_el") -> Path:
    suffix = "" if fan_layout == "az_el" else "_equal_area_guarded_D"
    return DATA_DIR / f"synthetic_truth_returns_{method}{suffix}_2-10MHz_100kHz.npz"


def chunk_path(method: str, start: int, stop: int, fan_layout: str = "az_el") -> Path:
    suffix = "" if fan_layout == "az_el" else "_equal_area_guarded_D"
    return DATA_DIR / f"synthetic_truth_{method}{suffix}_chunk_{start:02d}_{stop:02d}.npz"


def merge_chunks(method: str, anchor_stride: int, anchor_block_size: int,
                 anchor_elevation_stride: int,
                 anchor_optimizer: str, fan_layout: str,
                 outer_fraction: float, guard_seed_limit: int | None) -> None:
    chunk_ranges = [(start, min(start + 10, len(FREQUENCIES))) for start in range(0, len(FREQUENCIES), 10)]
    record_chunks = []
    counts = np.zeros((len(FREQUENCIES), 2), dtype=int)
    anchor_fan_directions = None
    full_fan_directions = None
    for start, stop in chunk_ranges:
        with np.load(chunk_path(method, start, stop, fan_layout), allow_pickle=False) as chunk:
            chunk_layout = str(chunk["vertical_fan_layout"]) if "vertical_fan_layout" in chunk.files else "az_el"
            if (str(chunk["method"]) != method or int(chunk["anchor_stride"]) != anchor_stride
                    or int(chunk["anchor_block_size"]) != anchor_block_size
                    or int(chunk["anchor_elevation_stride"]) != anchor_elevation_stride
                    or str(chunk["anchor_optimizer"]) != anchor_optimizer
                    or chunk_layout != fan_layout
                    or float(chunk["vertical_outer_ray_fraction"]) != outer_fraction
                    or int(chunk["vertical_guard_seed_limit"]) != (-1 if guard_seed_limit is None else guard_seed_limit)):
                raise ValueError(f"Inconsistent sweep settings in {chunk_path(method, start, stop, fan_layout)}")
            record_chunks.append(np.asarray(chunk["records"], dtype=float))
            counts += np.asarray(chunk["count_array"], dtype=int)
            chunk_fan_directions = int(chunk["anchor_fan_launch_directions"])
            chunk_full_directions = int(chunk["fan_launch_directions"])
            if anchor_fan_directions is not None and chunk_fan_directions != anchor_fan_directions:
                raise ValueError("Inconsistent anchor fan sizes across chunks")
            if full_fan_directions is not None and chunk_full_directions != full_fan_directions:
                raise ValueError("Inconsistent full fan sizes across chunks")
            anchor_fan_directions = chunk_fan_directions
            full_fan_directions = chunk_full_directions
    records = np.concatenate(record_chunks)
    output = output_path(method, fan_layout)
    np.savez_compressed(
        output,
        records=records,
        count_array=counts,
        frequencies_mhz=FREQUENCIES,
        method=np.array(method),
        vertical_fan_layout=np.array(fan_layout),
        vertical_outer_ray_fraction=np.array(outer_fraction),
        vertical_guard_seed_limit=np.array(-1 if guard_seed_limit is None else guard_seed_limit),
        anchor_stride=np.array(anchor_stride),
        anchor_block_size=np.array(anchor_block_size),
        anchor_elevation_stride=np.array(anchor_elevation_stride),
        anchor_optimizer=np.array(anchor_optimizer),
        fan_launch_directions=np.array(full_fan_directions),
        anchor_fan_launch_directions=np.array(anchor_fan_directions),
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
    parser.add_argument("--output", type=Path,
                        help="Output NPZ for a complete 2–10 MHz sweep")
    parser.add_argument("--method", choices=("adaptive", "dense"), default="adaptive")
    parser.add_argument("--density-scale", type=float, default=1.12,
                        help="Synthetic electron-density multiplier (default: 1.12)")
    parser.add_argument("--anchor-stride", type=int, default=5)
    parser.add_argument("--anchor-block-size", type=int, default=10)
    parser.add_argument("--anchor-elevation-stride", type=int, default=2)
    parser.add_argument("--anchor-optimizer", choices=("Powell", "Nelder-Mead"), default="Powell")
    parser.add_argument("--vertical-fan-layout", choices=("az_el", "equal_area_guarded"),
                        default="equal_area_guarded")
    parser.add_argument("--vertical-outer-ray-fraction", type=float,
                        help="Fraction of equal-area directions outside the near-nadir guard")
    parser.add_argument("--vertical-guard-seed-limit", type=int,
                        help="Additional near-nadir minima to optimize; 0 uses nearest-neighbor minima only, -1 uses all")
    args = parser.parse_args()
    if args.vertical_fan_layout == "equal_area_guarded":
        if args.vertical_outer_ray_fraction is None:
            args.vertical_outer_ray_fraction = 0.5
        if args.vertical_guard_seed_limit is None:
            args.vertical_guard_seed_limit = 4
    else:
        if args.vertical_outer_ray_fraction is None:
            args.vertical_outer_ray_fraction = 1.0
        if args.vertical_guard_seed_limit is None:
            args.vertical_guard_seed_limit = -1
    if args.anchor_stride < 1:
        parser.error("--anchor-stride must be positive")
    if args.anchor_block_size < 1:
        parser.error("--anchor-block-size must be positive")
    if args.anchor_elevation_stride < 1:
        parser.error("--anchor-elevation-stride must be positive")
    if not 0.0 < args.vertical_outer_ray_fraction <= 1.0:
        parser.error("--vertical-outer-ray-fraction must be in (0, 1]")
    if args.vertical_guard_seed_limit < -1:
        parser.error("--vertical-guard-seed-limit must be -1 or nonnegative")
    if not 0.0 < args.density_scale < 5.0:
        parser.error("--density-scale must be between 0 and 5")
    if args.density_scale != 1.12 and args.output is None:
        parser.error("non-default density scales require --output to keep each result separate")
    guard_seed_limit = None if args.vertical_guard_seed_limit == -1 else args.vertical_guard_seed_limit
    if args.vertical_fan_layout == "az_el":
        if args.vertical_outer_ray_fraction != 1.0 or guard_seed_limit is not None:
            parser.error("az_el fan does not accept guarded equal-area reductions")
    elif ((args.vertical_outer_ray_fraction, guard_seed_limit) != (0.5, 4)
          and args.output is None):
        parser.error("non-default guarded fans require --output to keep each result separate")
    if args.merge:
        if args.output is not None:
            parser.error("--output cannot be used with --merge")
        merge_chunks(args.method, args.anchor_stride, args.anchor_block_size,
                     args.anchor_elevation_stride,
                     args.anchor_optimizer, args.vertical_fan_layout,
                     args.vertical_outer_ray_fraction, guard_seed_limit)
        return
    if args.start is None or args.stop is None or not 0 <= args.start < args.stop <= len(FREQUENCIES):
        parser.error("choose a batch with --start and --stop between 0 and 81")
    if args.output is not None and (args.start != 0 or args.stop != FREQUENCIES.size):
        parser.error("--output requires the complete 0–81 frequency sweep")
    sweep_start = time.perf_counter()
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
            vertical_fan_layout=args.vertical_fan_layout,
            vertical_outer_ray_fraction=args.vertical_outer_ray_fraction,
            vertical_guard_seed_limit=guard_seed_limit,
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
        if args.method == "adaptive":
            unique_elevations = np.unique(case.fan_elevations_deg)
            retained = np.zeros(unique_elevations.size, dtype=bool)
            retained[:2] = True
            retained[2::args.anchor_elevation_stride] = True
            retained[-1] = True
            anchor_fan_directions = int(np.count_nonzero(np.isin(
                case.fan_elevations_deg, unique_elevations[retained])))
        else:
            anchor_fan_directions = int(case.fan_elevations_deg.size)
        grid = _apply_fit_params_to_grid(
            problem,
            problem.background_grids[0],
            IonosphereFitParams(
                density_scale=args.density_scale,
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
                    anchor_block_size=args.anchor_block_size,
                    anchor_elevation_stride=args.anchor_elevation_stride,
                    anchor_optimizer=args.anchor_optimizer,
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
        full_sweep = args.start == 0 and args.stop == FREQUENCIES.size
        destination = (args.output or output_path(args.method, args.vertical_fan_layout)) if full_sweep else chunk_path(
            args.method, args.start, args.stop, args.vertical_fan_layout)
        metadata = dict(
            records=records_array,
            count_array=counts,
            method=np.array(args.method),
            vertical_fan_layout=np.array(args.vertical_fan_layout),
            anchor_stride=np.array(args.anchor_stride),
            anchor_block_size=np.array(args.anchor_block_size),
            anchor_elevation_stride=np.array(args.anchor_elevation_stride),
            anchor_optimizer=np.array(args.anchor_optimizer),
            anchor_fan_launch_directions=np.array(anchor_fan_directions),
            fan_launch_directions=np.array(case.fan_elevations_deg.size),
            vertical_outer_ray_fraction=np.array(args.vertical_outer_ray_fraction),
            vertical_guard_seed_limit=np.array(args.vertical_guard_seed_limit),
            density_scale=np.array(args.density_scale),
            runtime_seconds=np.array(time.perf_counter() - sweep_start),
        )
        if full_sweep:
            metadata.update(
                frequencies_mhz=FREQUENCIES,
                homing_tolerance_m=np.array(config.homing_tolerance_m),
            )
        destination.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(destination, **metadata)
        print(destination, flush=True)


if __name__ == "__main__":
    main()
