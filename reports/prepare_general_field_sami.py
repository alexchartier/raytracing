"""Freeze seven SAMI soundings and synthetic in-situ measurements for inversion.

This is the only preparation step that reads the SAMI density truth. The
candidate builder and ionogram scorer consume only the two small JSON files,
the observed ionograms, and the independent PyIRI prior.
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
from python_raytrace.grid import load_ionosphere_grid_netcdf  # noqa: E402

CASE = ROOT / "reports/data/sami3_wave_20170111_0600"
OUTPUT = ROOT / "reports/data/general_field_sami3_wave"
INDICES = (1, 5, 9, 12, 14, 17, 20)


def main() -> None:
    os.umask(0o077)
    OUTPUT.mkdir(parents=True, exist_ok=True, mode=0o700)
    source = json.loads((CASE / "manifest.json").read_text())
    profiles = [row for row in source["profiles"] if row["index"] in INDICES]
    if tuple(row["index"] for row in profiles) != INDICES:
        raise ValueError("Missing saved truth ionogram position")
    grid = load_ionosphere_grid_netcdf(CASE / "truth_grid.nc")
    points = np.array([(row["latitude_deg"], row["longitude_deg"],
                        row["altitude_km"]) for row in profiles])
    measured = RegularGridInterpolator(
        (grid.latitudes_deg, grid.longitudes_deg, grid.altitudes_km),
        grid.iono_en_grid, bounds_error=True)(points)
    manifest = {key: value for key, value in source.items() if key != "profiles"}
    manifest["profiles"] = profiles
    manifest["profile_count"] = len(profiles)
    manifest["selection_observables"] = [
        "seven raw adaptive O/X ionograms", "seven local densities at 800 km"]
    observation = {
        "altitude_km": 800.0,
        "profile_indices": list(INDICES),
        "density_cm3_by_profile": measured.tolist(),
        "relative_uncertainty": 0.005,
        "provenance": "synthetic in-situ samples of SAMI truth; full density withheld from candidate selection",
    }
    for filename, data in (("manifest.json", manifest),
                           ("local_density_observations.json", observation)):
        path = OUTPUT / filename
        path.write_text(json.dumps(data, indent=2) + "\n")
        path.chmod(0o600)
    print(json.dumps({"indices": list(INDICES),
                      "local_density_cm3": measured.tolist()}, indent=2))


if __name__ == "__main__":
    main()
