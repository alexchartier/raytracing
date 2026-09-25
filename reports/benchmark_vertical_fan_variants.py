"""Check fast vertical fans on alternate synthetic ionospheres.

Use --full-sweep to save complete A/B/C/D accepted-return files for both
profiles. The default and --tune-guard retain the shorter benchmark windows.
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

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reports.benchmark_adaptive_homing import matched_returns  # noqa: E402
from reports.build_topside_validation_report import write_orbit_fixture  # noqa: E402
from python_raytrace.multisat_topside_inverse_demo import (  # noqa: E402
    IonosphereFitParams,
    TopsideInverseConfig,
    _apply_fit_params_to_grid,
    _equal_area_vertical_fan,
    _subset_problem_cases,
    build_inverse_problem,
    home_frequency_sweep_adaptive,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tune-guard", action="store_true",
                        help="Compare half-density fans with 1, 2, 4, or all extra guard minima")
    parser.add_argument("--full-sweep", action="store_true",
                        help="Run and save 2–10 MHz A/B/C/D ionograms for both profiles")
    args = parser.parse_args()
    if args.tune_guard and args.full_sweep:
        parser.error("--tune-guard and --full-sweep cannot be combined")
    windows = (("shape_bump", 2.0, 10.0), ("wave_3d", 2.0, 10.0)) if args.full_sweep else (
        ("shape_bump", 2.0, 2.4), ("wave_3d", 8.0, 8.4))
    variants = ((("A", None, None), ("B", 1.0, 0), ("C", 0.5, 0),
                 ("D", 0.5, 4)) if args.full_sweep else
                (("A", None, None), ("G1", 0.5, 1), ("G2", 0.5, 2),
                 ("G4", 0.5, 4), ("Gall", 0.5, None)) if args.tune_guard else
                (("A", None, None), ("B", 1.0, 0), ("C", 0.5, 0), ("Q", 0.25, 0)))
    output_name = ("vertical_fan_alternate_full.json" if args.full_sweep else
                   "vertical_fan_guard_tuning.json" if args.tune_guard else
                   "vertical_fan_cross_profile.json")
    output = ROOT / "reports" / "data" / output_name
    results = []
    with tempfile.TemporaryDirectory(prefix="vertical_fan_check_") as directory:
        temporary = Path(directory)
        orbit = temporary / "generated_orbit.nc"
        write_orbit_fixture(orbit)
        config = replace(
            TopsideInverseConfig(), ampere_file=orbit, planes=(1,),
            frequencies_mhz=(4.0, 5.0, 6.0), vertical_elevation_count=24,
            vertical_azimuth_step_deg=30.0, oblique_elevation_count=9,
            oblique_bearing_count=7, seed_max_candidates_per_frequency=32,
            homed_max_returns_per_frequency=64, grid_lat_step_deg=2.0,
            grid_lon_step_deg=4.0, grid_alt_step_km=20.0,
            d_region_model="none", grid_cache_path=temporary / "background.nc",
        )
        problem = _subset_problem_cases(build_inverse_problem(config), ["plane01_sv121_vertical"])
        case = problem.cases[0]
        background = problem.background_grids[0]
        for profile, start_mhz, stop_mhz in windows:
            if profile == "shape_bump":
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
            frequencies = np.round(np.arange(start_mhz, stop_mhz + 0.001, 0.1), 5)
            by_variant = {}
            for name, fraction, guard_limit in variants:
                if fraction is None:
                    variant_config = config
                    elevations, bearings = case.fan_elevations_deg, case.fan_bearings_deg
                else:
                    variant_config = replace(config, vertical_fan_layout="equal_area_guarded",
                                             vertical_outer_ray_fraction=fraction,
                                             vertical_guard_seed_limit=guard_limit)
                    elevations, bearings = _equal_area_vertical_fan(variant_config, guard_nadir_rows=3)
                by_cell = {}
                started = time.perf_counter()
                for mode in (1, -1):
                    return_sets = home_frequency_sweep_adaptive(
                        tx=case.tx_points[1], rx=case.rx_points[1], grid=grid,
                        fan_elevations_deg=elevations, fan_bearings_deg=bearings,
                        frequencies_mhz=frequencies, ox_mode=mode, config=variant_config,
                        range_min_km=150.0, anchor_stride=5, anchor_block_size=10,
                        anchor_elevation_stride=2,
                    )
                    by_cell.update({(i, mode): returns for i, returns in enumerate(return_sets)})
                by_variant[name] = by_cell
                runtime_seconds = time.perf_counter() - started
                reference = by_variant["A"]
                matched = 0
                missing = []
                for key, returns in by_cell.items():
                    count, absent = matched_returns(reference[key], returns)
                    matched += count
                    missing.extend({"frequency_mhz": float(frequencies[key[0]]),
                                    "mode": key[1], **value} for value in absent)
                entry = {
                    "profile": profile, "window_mhz": [start_mhz, stop_mhz],
                    "variant": name, "fan_directions": int(elevations.size),
                    "guard_seed_limit": guard_limit,
                    "runtime_seconds": runtime_seconds,
                    "reference_returns": sum(len(cell) for cell in reference.values()),
                    "variant_returns": sum(len(cell) for cell in by_cell.values()),
                    "matched_reference_returns": matched, "missing": missing,
                    "homing_tolerance_m": variant_config.homing_tolerance_m,
                }
                results.append(entry)
                print(entry, flush=True)
                if args.full_sweep:
                    unique_elevations = np.unique(elevations)
                    retained = np.zeros(unique_elevations.size, dtype=bool)
                    retained[:2] = True
                    retained[2::2] = True
                    retained[-1] = True
                    anchor_count = int(np.count_nonzero(np.isin(elevations, unique_elevations[retained])))
                    records = np.asarray([
                        (float(i), float(mode), ray.group_range_km, ray.miss_m, ray.absorption_db)
                        for mode in (1, -1) for i in range(frequencies.size)
                        for ray in by_cell[(i, mode)]
                    ], dtype=float).reshape(-1, 5)
                    destination = ROOT / "reports" / "data" / f"vertical_fan_{profile}_{name}_full.npz"
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    np.savez_compressed(
                        destination, records=records, frequencies_mhz=frequencies,
                        profile=np.array(profile), method=np.array("adaptive"),
                        vertical_fan_layout=np.array("az_el" if fraction is None else "equal_area_guarded"),
                        vertical_outer_ray_fraction=np.array(1.0 if fraction is None else fraction),
                        vertical_guard_seed_limit=np.array(-1 if guard_limit is None else guard_limit),
                        anchor_fan_launch_directions=np.array(anchor_count),
                        fan_launch_directions=np.array(elevations.size),
                        homing_tolerance_m=np.array(variant_config.homing_tolerance_m),
                        runtime_seconds=np.array(runtime_seconds),
                    )
                    print(destination, flush=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(results, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
