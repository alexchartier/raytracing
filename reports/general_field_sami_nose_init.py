"""Make a peak-frequency-initialized SAMI candidate without density truth.

Relative O/X nose shifts against the independent PyIRI forward ionograms give
the required log peak-density changes. A smooth along-track basis interpolates
those seven observations; the local 800 km measurements anchor the upper tail.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from ionogram_metrics import Ionogram  # noqa: E402
from python_raytrace.general_field_inverse import (  # noqa: E402
    LocalDensity, TrackLogDensityBasis, candidate_grid, sampled_local_density,
)
from python_raytrace.grid import (  # noqa: E402
    load_ionosphere_grid_netcdf, save_ionosphere_grid_netcdf,
)

CASE = ROOT / "reports/data/sami3_wave_20170111_0600"
RUN = ROOT / "reports/data/general_field_sami3_wave"


def main() -> None:
    os.umask(0o077)
    plan = json.loads((RUN / "plan.json").read_text())
    profiles = json.loads((ROOT / plan["manifest"]).read_text())["profiles"]
    indices = [row["index"] for row in profiles]
    locations = np.array([(row["latitude_deg"], row["longitude_deg"])
                          for row in profiles])
    observation = json.loads((ROOT / plan["local_density_observations"]).read_text())
    local = LocalDensity(locations, observation["altitude_km"],
                         np.array(observation["density_cm3_by_profile"]),
                         observation["relative_uncertainty"])
    basis = TrackLogDensityBasis.from_locations(
        locations, spatial_centers=plan["spatial_centers"],
        altitude_offsets_km=tuple(plan["vertical_offsets_km"]),
        altitude_widths_km=tuple(plan["vertical_widths_km"]))
    prior = load_ionosphere_grid_netcdf(ROOT / plan["prior_grid"])
    target = []
    raw = []
    for index in indices:
        observed = Ionogram.read(CASE / "vertical_truth" / f"ionogram_{index:02d}.npz")
        predicted = Ionogram.read(CASE / "vertical_prior" / f"ionogram_{index:02d}.npz")
        pairs = [(observed.nose(mode), predicted.nose(mode)) for mode in (1, -1)]
        shifts = [2 * np.log(a / b) for a, b in pairs
                  if a is not None and b is not None and a < 9.9 and b < 9.9]
        if not shifts:
            raise ValueError(f"No usable O/X nose for {index}")
        target.append(float(np.mean(shifts)))
        raw.append({"index": index, "observed_and_prior_noses_mhz": pairs,
                    "target_log_peak_density_change": target[-1]})
    distance = basis.frame.sample_distance_km
    centers = basis.center_distance_km
    spatial = np.exp(-.5 * ((distance[:, None] - centers) / basis.spatial_width_km) ** 2)
    spatial /= spatial.sum(axis=1, keepdims=True)
    second = np.diff(np.eye(len(centers)), n=2, axis=0)
    ridge = .008
    system = spatial.T @ spatial + ridge * (second.T @ second + .2 * np.eye(len(centers)))
    fitted = np.linalg.solve(system, spatial.T @ np.asarray(target))
    coefficients = np.zeros(basis.coefficient_shape)
    peak_index = int(np.argmin(abs(basis.altitude_offsets_km)))
    coefficients[:, peak_index] = fitted
    candidate = candidate_grid(prior, basis, coefficients, local_density=local,
                               maximum_log_correction=plan["maximum_log_correction"])
    output = RUN / "nose_init"
    output.mkdir(mode=0o700, parents=True, exist_ok=True)
    grid_path = output / "grid.nc"
    save_ionosphere_grid_netcdf(grid_path, candidate)
    grid_path.chmod(0o600)
    diagnostic = {
        "selection_uses_truth_density": False,
        "source_ionograms": ["SAMI O/X observed", "PyIRI O/X prior"],
        "source_local_density": "seven synthetic 800 km in-situ samples",
        "nose_observations": raw,
        "basis_coefficients": coefficients.tolist(),
        "fitted_log_peak_density_change": (spatial @ fitted).tolist(),
        "maximum_local_density_error_percent": float(np.max(abs(
            100 * (sampled_local_density(candidate, local)
                   / local.electron_density_cm3 - 1)))),
    }
    derivation = output / "derivation.json"
    derivation.write_text(json.dumps(diagnostic, indent=2) + "\n")
    derivation.chmod(0o600)
    plan["candidates"] = [row for row in plan["candidates"]
                          if row["name"] != "nose_init"]
    plan["candidates"].append({"name": "nose_init",
                               "grid": str(grid_path.relative_to(ROOT)),
                               "coefficients": coefficients.tolist(),
                               "maximum_local_density_error_percent":
                                   diagnostic["maximum_local_density_error_percent"]})
    (RUN / "plan.json").write_text(json.dumps(plan, indent=2) + "\n")
    (RUN / "plan.json").chmod(0o600)
    print(json.dumps({"target": target,
                      "fitted": diagnostic["fitted_log_peak_density_change"]}, indent=2))


if __name__ == "__main__":
    main()
