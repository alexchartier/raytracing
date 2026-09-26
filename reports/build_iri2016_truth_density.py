"""Build an independent IRI-2016 density grid for the D sounder geometry.

Requires the public ``iri2016`` Python package and its Fortran compiler.
Only electron density comes from IRI-2016; the D ray tracer supplies the
same geometry and geomagnetic field used for the candidate ionograms.
"""

from __future__ import annotations

import argparse
import datetime as dt
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

import numpy as np
from iri2016 import IRI

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reports.build_topside_validation_report import write_orbit_fixture  # noqa: E402
from python_raytrace.multisat_topside_inverse_demo import (  # noqa: E402
    TopsideInverseConfig, _subset_problem_cases, build_inverse_problem,
)


def _repair_density(altitudes_km: np.ndarray, density_m3: np.ndarray) -> np.ndarray:
    valid = np.isfinite(density_m3) & (density_m3 > 0.0)
    if np.count_nonzero(valid) < 2:
        raise ValueError("IRI-2016 supplied fewer than two valid density samples")
    result = density_m3.copy()
    invalid = ~valid
    if np.any(invalid):
        log_density = np.log(density_m3[valid])
        result[invalid] = np.exp(np.interp(altitudes_km[invalid], altitudes_km[valid],
                                           log_density))
        first = np.flatnonzero(valid)[0]
        if first > 0:
            second = np.flatnonzero(valid)[1]
            slope = (np.log(density_m3[second]) - np.log(density_m3[first])) / (
                altitudes_km[second] - altitudes_km[first])
            result[:first] = np.exp(np.log(density_m3[first])
                                    + slope * (altitudes_km[:first] - altitudes_km[first]))
    return np.maximum(result, 1e-6)


def _local_axes() -> tuple[np.ndarray, np.ndarray, np.ndarray, dt.datetime]:
    with tempfile.TemporaryDirectory(prefix="iri2016_d_truth_") as directory:
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
        background = problem.background_grids[0]
        return (np.asarray(background.latitudes_deg, dtype=float),
                np.asarray(background.longitudes_deg, dtype=float),
                np.asarray(background.altitudes_km, dtype=float), case.times_utc[1])


def build(output: Path, axes_path: Path | None = None) -> None:
    if axes_path is None:
        lats, lons, alts, when = _local_axes()
    else:
        with np.load(axes_path, allow_pickle=False) as axes:
            lats = np.asarray(axes["latitudes_deg"], dtype=float)
            lons = np.asarray(axes["longitudes_deg"], dtype=float)
            alts = np.asarray(axes["altitudes_km"], dtype=float)
            when = dt.datetime.fromisoformat(str(axes["time_utc"]))
    step = float(alts[1] - alts[0])
    densities_cm3 = np.empty((len(lats), len(lons), len(alts)), dtype=float)
    invalid_count = 0
    f107 = []
    for lat_index, lat in enumerate(lats):
        for lon_index, lon in enumerate(lons):
            profile = IRI(when, [float(alts[0]), float(alts[-1]), step],
                          float(lat), float(lon))
            if not np.allclose(profile.alt_km.values, alts):
                raise ValueError("IRI-2016 altitude axis differs from the D grid")
            density_m3 = np.asarray(profile["ne"].values, dtype=float)
            invalid_count += int(np.count_nonzero(~np.isfinite(density_m3) | (density_m3 <= 0)))
            densities_cm3[lat_index, lon_index] = _repair_density(alts, density_m3) / 1e6
            f107.append(float(profile.attrs["f107"]))
        print(f"{lat_index + 1}/{len(lats)} latitude rows", flush=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, latitudes_deg=lats, longitudes_deg=lons,
                        altitudes_km=alts, electron_density_cm3=densities_cm3,
                        model=np.array("IRI-2016 1.11.1"),
                        time_utc=np.array(when.isoformat()),
                        invalid_samples_repaired=np.array(invalid_count),
                        f107_min=np.array(min(f107)), f107_max=np.array(max(f107)))
    print(output)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--axes", type=Path,
                        help="D case axes exported by the same runtime that will ray trace")
    args = parser.parse_args()
    build(args.output, args.axes)


if __name__ == "__main__":
    main()
