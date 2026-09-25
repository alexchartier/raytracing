"""Compare adaptive and full-fan homing on distinct synthetic ionospheres.

Run one profile at a time, for example:
``python3 reports/benchmark_adaptive_homing.py --profile scale_09``.
The result is written to reports/data/homing_benchmark_<profile>.json.
"""

from __future__ import annotations

import argparse
import json
import sys
import tempfile
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

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

WINDOWS_MHZ = ((2.0, 2.4), (4.0, 4.4), (8.0, 8.4))


def matched_returns(dense, adaptive, *, range_tolerance_km: float = 1.0,
                    angle_tolerance_deg: float = 0.5):
    dense_ranges = np.asarray([return_.group_range_km for return_ in dense], dtype=float)
    adaptive_ranges = np.asarray([return_.group_range_km for return_ in adaptive], dtype=float)
    def launch_directions(returns):
        elevation = np.deg2rad([float(return_.ray.path["initial_elev"]) for return_ in returns])
        bearing = np.deg2rad([float(return_.ray.path["initial_bearing"]) for return_ in returns])
        return np.column_stack((np.cos(elevation) * np.sin(bearing),
                                np.cos(elevation) * np.cos(bearing), np.sin(elevation)))

    def missing(indices):
        return [{"range_km": float(dense_ranges[index]),
                 "elevation_deg": float(dense[index].ray.path["initial_elev"]),
                 "bearing_deg": float(dense[index].ray.path["initial_bearing"])}
                for index in indices]

    if not dense_ranges.size or not adaptive_ranges.size:
        return 0, missing(range(len(dense)))
    distances = np.abs(dense_ranges[:, None] - adaptive_ranges[None, :])
    directions_dense = launch_directions(dense)
    directions_adaptive = launch_directions(adaptive)
    angles = np.rad2deg(np.arccos(np.clip(directions_dense @ directions_adaptive.T, -1.0, 1.0)))
    valid = (distances <= range_tolerance_km) & (angles <= angle_tolerance_deg)
    cost = np.where(valid, distances / range_tolerance_km + angles / angle_tolerance_deg, 1e6)
    row_indices, column_indices = linear_sum_assignment(cost)
    matched_rows = {
        int(row) for row, column in zip(row_indices, column_indices)
        if valid[row, column]
    }
    return len(matched_rows), missing(index for index in range(len(dense)) if index not in matched_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", choices=("scale_09", "shape_bump", "wave_3d"), required=True)
    args = parser.parse_args()
    output = ROOT / "reports" / "data" / f"homing_benchmark_{args.profile}.json"

    with tempfile.TemporaryDirectory(prefix="topside_homing_benchmark_") as directory:
        temporary = Path(directory)
        orbit = temporary / "generated_orbit.nc"
        write_orbit_fixture(orbit)
        config = replace(
            TopsideInverseConfig(),
            ampere_file=orbit,
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
        background = problem.background_grids[0]
        if args.profile == "scale_09":
            grid = _apply_fit_params_to_grid(problem, background, IonosphereFitParams(
                density_scale=0.9, hmf2_shift_km=0.0, wave_amplitude_fraction=0.0,
                wave_phase_rad=0.0, wave_bearing_deg=0.0,
            ))
        elif args.profile == "shape_bump":
            heights = np.asarray(background.altitudes_km, dtype=float)
            factor = 1.10 + 0.10 * np.exp(-0.5 * ((heights - 320.0) / 55.0) ** 2)
            grid = replace(
                background,
                iono_en_grid=np.asarray(background.iono_en_grid) * factor[None, None, :],
                iono_en_grid_5=np.asarray(background.iono_en_grid_5) * factor[None, None, :],
            )
        else:
            grid = _apply_fit_params_to_grid(problem, background, IonosphereFitParams(
                density_scale=1.0, hmf2_shift_km=10.0, wave_amplitude_fraction=0.18,
                wave_phase_rad=0.6, wave_bearing_deg=45.0,
            ))

        result = {
            "profile": args.profile,
            "fan_directions": int(case.fan_elevations_deg.size),
            "homing_tolerance_m": float(config.homing_tolerance_m),
            "anchor_stride": 5,
            "range_match_tolerance_km": 1.0,
            "angle_match_tolerance_deg": 0.5,
            "windows": [],
        }
        for first_mhz, last_mhz in WINDOWS_MHZ:
            frequencies = np.round(np.arange(first_mhz, last_mhz + 0.001, 0.1), 5)
            dense_by_cell = {}
            adaptive_by_cell = {}
            dense_seconds = 0.0
            adaptive_seconds = 0.0
            for mode in (1, -1):
                start = time.perf_counter()
                for frequency_index, frequency_mhz in enumerate(frequencies):
                    dense_by_cell[(frequency_index, mode)] = tuple(
                        solution for solution in _home_frequency_returns(
                            PointToPointRayTracer(),
                            tx=case.tx_points[1], rx=case.rx_points[1], grid=grid,
                            fan_elevations_deg=case.fan_elevations_deg,
                            fan_bearings_deg=case.fan_bearings_deg,
                            frequency_mhz=float(frequency_mhz), ox_mode=mode,
                            config=config,
                        ) if solution.group_range_km >= 150.0
                    )
                dense_seconds += time.perf_counter() - start

                start = time.perf_counter()
                adaptive_sets = home_frequency_sweep_adaptive(
                    tx=case.tx_points[1], rx=case.rx_points[1], grid=grid,
                    fan_elevations_deg=case.fan_elevations_deg,
                    fan_bearings_deg=case.fan_bearings_deg,
                    frequencies_mhz=frequencies, ox_mode=mode, config=config,
                    range_min_km=150.0, anchor_stride=5,
                )
                adaptive_seconds += time.perf_counter() - start
                for frequency_index, solutions in enumerate(adaptive_sets):
                    adaptive_by_cell[(frequency_index, mode)] = solutions

            dense_count = sum(len(cell) for cell in dense_by_cell.values())
            adaptive_count = sum(len(cell) for cell in adaptive_by_cell.values())
            matched_count = 0
            misses = []
            for frequency_index, frequency_mhz in enumerate(frequencies):
                for mode in (1, -1):
                    matched, missing_ranges = matched_returns(
                        dense_by_cell[(frequency_index, mode)],
                        adaptive_by_cell[(frequency_index, mode)],
                    )
                    matched_count += matched
                    misses.extend({"frequency_mhz": float(frequency_mhz), "mode": mode, **value}
                                  for value in missing_ranges)
            window_result = {
                "first_mhz": first_mhz,
                "last_mhz": last_mhz,
                "dense_seconds": dense_seconds,
                "adaptive_seconds": adaptive_seconds,
                "speedup": dense_seconds / adaptive_seconds,
                "dense_returns": dense_count,
                "adaptive_returns": adaptive_count,
                "matched_dense_returns": matched_count,
                "misses": misses,
            }
            result["windows"].append(window_result)
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(json.dumps(result, indent=2) + "\n")
            print(window_result, flush=True)
    print(output)


if __name__ == "__main__":
    main()
