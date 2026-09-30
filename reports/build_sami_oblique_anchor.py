"""Extend the varying-topside fit to all oblique spacecraft in-situ positions."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from python_raytrace.general_field_inverse import LocalDensity  # noqa: E402
from python_raytrace.grid import load_ionosphere_grid_netcdf, save_ionosphere_grid_netcdf  # noqa: E402
from python_raytrace.monotone_topside import monotone_topside_grid  # noqa: E402

CASE = ROOT / "reports/data/sami3_wave_20170111_0600"
VERTICAL = ROOT / "reports/data/sami3_monotone_topside"
RUN = ROOT / "reports/data/sami3_monotone_oblique"


def main() -> None:
    os.umask(0o077)
    insitu = json.loads((RUN / "insitu_1km.json").read_text())
    local = LocalDensity(np.asarray(insitu["lat_lon_deg"]),
                         insitu["altitude_km"],
                         np.asarray(insitu["density_cm3"]),
                         insitu["relative_uncertainty"])
    plan = json.loads((VERTICAL / "plan.json").read_text())
    selected = json.loads((VERTICAL / "selection.json").read_text())
    if selected["chosen_candidate"] != "varying_g04":
        raise ValueError("Expected frozen varying_g04 vertical starting field")
    old = next(row for row in plan["candidates"]
               if row["name"] == selected["chosen_candidate"])
    profiles = json.loads((CASE / "manifest.json").read_text())["profiles"]
    indices = plan["observed_profile_indices"]
    sounding_lat = np.array([profiles[i - 1]["latitude_deg"] for i in indices])
    grid = load_ionosphere_grid_netcdf(ROOT / old["grid"])
    candidate = monotone_topside_grid(
        grid, local, sounding_lat,
        np.asarray(old["fraction_at_quarter_by_profile"]),
        np.asarray(old["fraction_at_three_fifths_by_profile"]),
        np.zeros(len(indices)),
    )
    target = RUN / "extended_anchor"
    target.mkdir(mode=0o700, parents=True, exist_ok=True)
    output = target / "grid.nc"
    save_ionosphere_grid_netcdf(output, candidate)
    output.chmod(0o600)
    result = {
        "source_grid": old["grid"],
        "candidate_grid": str(output.relative_to(ROOT)),
        "dense_insitu": str((RUN / "insitu_1km.json").relative_to(ROOT)),
        "source_vertical_selection": str((VERTICAL / "selection.json").relative_to(ROOT)),
        "oblique_truth_ionograms": str((CASE / "oblique_wave_600km/truth").relative_to(ROOT)),
        "oblique_geometry": "600 km along-track, both spacecraft at 800 km",
        "observed_indices": indices,
        "density_truth_used_in_selection": False,
    }
    path = RUN / "plan.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    path.chmod(0o600)
    print(json.dumps({"grid": str(output), "insitu_samples": len(insitu["density_cm3"])},
                     indent=2))


if __name__ == "__main__":
    main()
