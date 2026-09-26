"""Verify the independent IRI-2016 density input and blinded ionogram.

Requires the public ``iri2016`` package used to create the grid. This audit
compares a fresh direct IRI call with the saved grid above the sounder.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
from pathlib import Path

import numpy as np
from iri2016 import IRI
from scipy.interpolate import RegularGridInterpolator


def audit(truth_ionogram: Path, truth_grid: Path, background_grid: Path,
          candidate_ionogram: Path, output: Path) -> None:
    with np.load(truth_ionogram, allow_pickle=False) as observed:
        truth_source = str(observed["density_source"])
        forbidden = ("density_scale", "hmf2_shift_km", "f2_width_scale",
                     "wave_amplitude_fraction", "wave_phase_rad", "wave_bearing_deg")
        exposed = [key for key in forbidden if key in observed.files]
    with np.load(candidate_ionogram, allow_pickle=False) as candidate:
        candidate_source = str(candidate["density_source"])
    with np.load(truth_grid, allow_pickle=False) as grid, np.load(background_grid, allow_pickle=False) as background:
        for axis in ("latitudes_deg", "longitudes_deg", "altitudes_km"):
            if not np.array_equal(grid[axis], background[axis]):
                raise ValueError(f"The truth and candidate {axis} axes differ")
        latitude = float(background["tx_lat_deg"])
        longitude = float(background["tx_lon_deg"])
        altitudes = np.asarray(grid["altitudes_km"], dtype=float)
        interpolated = RegularGridInterpolator(
            (grid["latitudes_deg"], grid["longitudes_deg"]),
            grid["electron_density_cm3"], bounds_error=True,
        )([[latitude, longitude]])[0]
        when = dt.datetime.fromisoformat(str(grid["time_utc"]))
        direct = np.asarray(IRI(when, [float(altitudes[0]), float(altitudes[-1]),
                                      float(altitudes[1] - altitudes[0])],
                                latitude, longitude)["ne"].values, dtype=float) / 1e6
        repaired_samples = int(grid["invalid_samples_repaired"])
        grid_samples = int(grid["electron_density_cm3"].size)
    if truth_source != "IRI-2016 1.11.1" or candidate_source != "PyIRI" or exposed:
        raise ValueError("The truth source, candidate source, or blinding check failed")
    mask = (altitudes >= 150) & (altitudes <= 600)
    if not np.all(np.isfinite(direct[mask])):
        raise ValueError("Direct IRI profile has invalid samples in the comparison band")
    discrepancy = float(np.sqrt(np.mean((interpolated[mask] - direct[mask]) ** 2))
                        / np.max(direct[mask]))
    result = {
        "truth_density_source": truth_source,
        "candidate_density_source": candidate_source,
        "truth_parameter_fields_exposed_to_fit": exposed,
        "truth_candidate_grid_axes_identical": True,
        "direct_iri_profile_vs_interpolated_grid_rms_relative_to_peak_150_to_600_km": discrepancy,
        "direct_iri_invalid_samples_150_to_600_km": int(np.count_nonzero(~np.isfinite(direct[mask]))),
        "truth_grid_repaired_samples": repaired_samples,
        "truth_grid_total_samples": grid_samples,
        "time_utc": when.isoformat(),
        "sounder_latitude_deg": latitude,
        "sounder_longitude_deg": longitude,
        "shared_components": ["geometry", "geomagnetic field", "ray tracer", "homing algorithm"],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("truth_ionogram", "truth_grid", "background_grid", "candidate_ionogram", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    audit(args.truth_ionogram, args.truth_grid, args.background_grid,
          args.candidate_ionogram, args.output)


if __name__ == "__main__":
    main()
