"""Export the D sounder grid axes from the runtime that will ray trace."""

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
    TopsideInverseConfig, _subset_problem_cases, build_inverse_problem,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="d_case_axes_") as directory:
        temp = Path(directory)
        orbit = temp / "orbit.nc"
        write_orbit_fixture(orbit)
        config = replace(
            TopsideInverseConfig(), ampere_file=orbit, planes=(1,),
            frequencies_mhz=(4.0, 5.0, 6.0), vertical_elevation_count=24,
            vertical_fan_layout="equal_area_guarded", vertical_outer_ray_fraction=0.5,
            vertical_guard_seed_limit=4, vertical_azimuth_step_deg=30.0,
            oblique_elevation_count=9, oblique_bearing_count=7,
            seed_max_candidates_per_frequency=32, homed_max_returns_per_frequency=64,
            grid_lat_step_deg=2.0, grid_lon_step_deg=4.0, grid_alt_step_km=20.0,
            d_region_model="none", grid_cache_path=temp / "background.nc",
        )
        problem = _subset_problem_cases(build_inverse_problem(config), ["plane01_sv121_vertical"])
        case = problem.cases[0]
        grid = problem.background_grids[0]
        args.output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.output,
                            latitudes_deg=np.asarray(grid.latitudes_deg, dtype=float),
                            longitudes_deg=np.asarray(grid.longitudes_deg, dtype=float),
                            altitudes_km=np.asarray(grid.altitudes_km, dtype=float),
                            pyiri_density_cm3=np.asarray(grid.iono_en_grid, dtype=float),
                            time_utc=np.array(case.times_utc[1].isoformat()),
                            tx_lat_deg=np.array(case.tx_points[1].lat_deg),
                            tx_lon_deg=np.array(case.tx_points[1].lon_deg))
        print(args.output)


if __name__ == "__main__":
    main()
