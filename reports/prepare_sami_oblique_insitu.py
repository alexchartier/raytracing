"""Save synthetic 800 km in-situ samples along both oblique spacecraft tracks."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from python_raytrace.general_field_inverse import EARTH_RADIUS_KM  # noqa: E402
from python_raytrace.grid import load_ionosphere_grid_netcdf  # noqa: E402

CASE = ROOT / "reports/data/sami3_wave_20170111_0600"
RUN = ROOT / "reports/data/sami3_monotone_oblique"


def main() -> None:
    os.umask(0o077)
    RUN.mkdir(parents=True, exist_ok=True, mode=0o700)
    os.environ["SOUNDER_CASE_ROOT"] = str(CASE)
    from oblique_wave_pass import endpoints  # noqa: E402

    first_tx, _ = endpoints(1)
    _, last_rx = endpoints(20)
    span_km = EARTH_RADIUS_KM * np.deg2rad(last_rx.lat_deg - first_tx.lat_deg)
    distance = np.r_[np.arange(0.0, span_km, 1.0), span_km]
    latitude = first_tx.lat_deg + np.rad2deg(distance / EARTH_RADIUS_KM)
    locations = np.column_stack((latitude, np.full(len(latitude), first_tx.lon_deg)))
    grid = load_ionosphere_grid_netcdf(CASE / "truth_grid.nc")
    density = RegularGridInterpolator(
        (grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km),
        grid.iono_en_grid, bounds_error=True)(
            np.column_stack((locations, np.full(len(latitude), first_tx.alt_km))))
    result = {
        "provenance": "synthetic 800 km density sampled from SAMI along both satellite tracks",
        "altitude_km": first_tx.alt_km,
        "spacing_km": 1.0,
        "distance_km": distance.tolist(),
        "lat_lon_deg": locations.tolist(),
        "density_cm3": density.tolist(),
        "relative_uncertainty": 0.005,
        "endpoints_latitude_deg": [first_tx.lat_deg, last_rx.lat_deg],
    }
    path = RUN / "insitu_1km.json"
    path.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    path.chmod(0o600)
    print(json.dumps({"samples": len(distance),
                      "maximum_spacing_km": float(np.max(np.diff(distance))),
                      "latitudes_deg": result["endpoints_latitude_deg"]}, indent=2))


if __name__ == "__main__":
    main()
