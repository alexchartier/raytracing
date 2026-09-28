"""Build shared topside-width candidates from the ionogram-selected wave grid.

The pilot grids perturb width uniformly. A later fit uses their traced response
to estimate a smooth latitude update; truth density is never read here.
"""

from __future__ import annotations

import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from python_raytrace.absorption import effective_collision_frequency  # noqa: E402
from python_raytrace.grid import load_ionosphere_grid_netcdf, save_ionosphere_grid_netcdf  # noqa: E402

DATA = ROOT / "reports/data"
START_GRID = DATA / "lat_wave_joint_round2/peak_full_grid.nc"
OUT = DATA / "lat_wave_width_round3"


def save_candidate(name: str, fractions: np.ndarray) -> dict:
    grid = load_ionosphere_grid_netcdf(START_GRID)
    z = np.asarray(grid.altitudes_km, dtype=float)
    density = np.asarray(grid.iono_en_grid, dtype=float)
    if fractions.shape != (len(grid.latitudes_deg),):
        raise ValueError("Expected one width fraction per latitude grid row")
    if np.max(abs(fractions)) > .15:
        raise ValueError("Topside width update exceeds 15%")
    modified = np.empty_like(density)
    for i, fraction in enumerate(fractions):
        for j in range(len(grid.longitudes_deg)):
            profile = density[i, j]
            peak_z = z[np.argmax(profile)]
            source_z = np.where(z > peak_z,
                                peak_z + (z - peak_z) / (1.0 + fraction), z)
            modified[i, j] = np.interp(source_z, z, profile,
                                       left=profile[0], right=profile[-1])
    if np.any(modified <= 0) or not np.all(np.isfinite(modified)):
        raise ValueError("Invalid candidate density")
    collision = grid.collision_freq
    if (grid.electron_temp_k is not None and grid.ion_temp_k is not None
            and grid.neutral_species_cm3 is not None):
        collision = effective_collision_frequency(
            grid.electron_temp_k, grid.ion_temp_k,
            modified * 1e6, grid.neutral_species_cm3)
    candidate_grid = replace(grid, iono_en_grid=modified,
                             iono_en_grid_5=modified, collision_freq=collision)
    OUT.mkdir(parents=True, exist_ok=True)
    grid_path = OUT / f"{name}_grid.nc"
    density_path = OUT / f"{name}_density.npz"
    save_ionosphere_grid_netcdf(grid_path, candidate_grid)
    np.savez_compressed(density_path, latitudes_deg=grid.latitudes_deg,
                        longitudes_deg=grid.longitudes_deg, altitudes_km=z,
                        electron_density_cm3=modified,
                        model=np.array(f"Joint wave retrieval, topside width {name}"))
    return {"name": name, "grid": str(grid_path.relative_to(ROOT)),
            "density": str(density_path.relative_to(ROOT)),
            "width_fraction_min": float(np.min(fractions)),
            "width_fraction_max": float(np.max(fractions))}


def build_pilots() -> None:
    grid = load_ionosphere_grid_netcdf(START_GRID)
    n = len(grid.latitudes_deg)
    candidates = [save_candidate("width_plus5", np.full(n, .05)),
                  save_candidate("width_minus5", np.full(n, -.05))]
    result = {"selection_uses_truth_density": False,
              "starting_grid": str(START_GRID.relative_to(ROOT)),
              "description": "Uniform +/-5% topside altitude-width pilots",
              "candidates": candidates}
    (OUT / "pilots.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    build_pilots()
