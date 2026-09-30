"""Save synthetic 800 km in-situ density every kilometre on the SAMI pass.

This preparation is the only new step that reads the full SAMI density grid.
The retrieval consumes the resulting one-dimensional observations, never the
three-dimensional truth field.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from python_raytrace.general_field_inverse import EARTH_RADIUS_KM  # noqa: E402
from python_raytrace.grid import load_ionosphere_grid_netcdf  # noqa: E402

CASE = ROOT / "reports/data/sami3_wave_20170111_0600"
RUN = ROOT / "reports/data/sami3_monotone_topside"


def main() -> None:
    os.umask(0o077)
    RUN.mkdir(parents=True, exist_ok=True, mode=0o700)
    grid = load_ionosphere_grid_netcdf(CASE / "truth_grid.nc")
    span_km = EARTH_RADIUS_KM * np.deg2rad(20.0)
    distance = np.r_[np.arange(0.0, span_km, 1.0), span_km]
    latitude = -70.0 + np.rad2deg(distance / EARTH_RADIUS_KM)
    locations = np.column_stack((latitude, np.full(len(latitude), 190.0)))
    points = np.column_stack((locations, np.full(len(latitude), 800.0)))
    observed = RegularGridInterpolator(
        (grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km),
        grid.iono_en_grid, bounds_error=True)(points)
    result = {
        "provenance": "synthetic 800 km in-situ track sampled from remapped SAMI density",
        "altitude_km": 800.0,
        "spacing_km": 1.0,
        "distance_km": distance.tolist(),
        "lat_lon_deg": locations.tolist(),
        "density_cm3": observed.tolist(),
        "relative_uncertainty": 0.005,
    }
    path = RUN / "insitu_1km.json"
    path.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    path.chmod(0o600)
    print(json.dumps({"samples": len(distance),
                      "maximum_spacing_km": float(np.max(np.diff(distance))),
                      "density_range_cm3": [float(observed.min()),
                                           float(observed.max())]}, indent=2))


if __name__ == "__main__":
    main()
