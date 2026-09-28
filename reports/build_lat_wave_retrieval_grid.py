"""Build the ionogram-selected PyIRI candidate grid for the latitude pass."""

from __future__ import annotations

import datetime as dt
import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from python_raytrace.absorption import effective_collision_frequency  # noqa: E402
from python_raytrace.grid import build_pyiri_grid_from_axes, save_ionosphere_grid_netcdf  # noqa: E402

DATA = ROOT / "reports/data"
PRIOR = DATA / "d_inverse_iri_peak_summary.json"
BACKGROUND = DATA / "d_inverse_iri_pyiri_background.npz"
INITIAL_FIT = DATA / "lat_wave_pass_initial_fit.json"
DENSITY_OUTPUT = DATA / "lat_wave_pass_retrieved_density.npz"
GRID_OUTPUT = DATA / "lat_wave_pass_retrieved_forward_grid.nc"


def build(fit_path: Path = INITIAL_FIT, density_output: Path = DENSITY_OUTPUT,
          grid_output: Path = GRID_OUTPUT) -> None:
    fit = json.loads(fit_path.read_text())
    prior = json.loads(PRIOR.read_text())["retrieved"]
    with np.load(BACKGROUND, allow_pickle=False) as source:
        latitudes = np.asarray(source["latitudes_deg"], dtype=float)
        longitudes = np.asarray(source["longitudes_deg"], dtype=float)
        altitudes = np.asarray(source["altitudes_km"], dtype=float)
        pyiri = np.asarray(source["pyiri_density_cm3"], dtype=float)
        when = dt.datetime.fromisoformat(str(source["time_utc"]))
    width = float(prior["f2_width_scale"])
    top_ratio = float(prior["topside_width_ratio"])
    shift = float(prior["hmf2_shift_km"])
    peak_fraction = float(prior["peak_perturbation_fraction"])
    peak_width = float(prior["peak_width_km"])
    density = np.empty_like(pyiri)
    for i in range(len(latitudes)):
        for j in range(len(longitudes)):
            source = pyiri[i, j]
            peak = float(altitudes[np.argmax(source)])
            side_width = np.where(altitudes >= peak, width * top_ratio, width)
            mapped = peak + (altitudes - peak) / side_width
            stretched = np.interp(mapped, altitudes, source,
                                  left=source[0], right=source[-1])
            shifted = np.interp(altitudes - shift, altitudes, stretched,
                                left=stretched[0], right=stretched[-1])
            if peak_fraction:
                shifted_peak = float(altitudes[np.argmax(shifted)])
                shifted *= 1.0 + peak_fraction * np.exp(
                    -0.5 * ((altitudes - shifted_peak) / peak_width) ** 2)
            density[i, j] = float(prior["density_scale"]) * shifted
    center_lat = float(np.mean(fit["pass_latitudes_deg"]))
    x_km = 6371.0088 * np.deg2rad(latitudes - center_lat)
    x_norm = x_km / max(abs(6371.0088 * np.deg2rad(
        np.asarray(fit["pass_latitudes_deg"]) - center_lat)))
    local_peak_altitudes = altitudes[np.argmax(density, axis=2)]
    envelope = np.exp(-0.5 * ((altitudes[None, None, :] - local_peak_altitudes[:, :, None])
                              / float(fit["wave_vertical_sigma_km_assumed"])) ** 2)
    wave = np.cos(2 * np.pi * x_km / float(fit["wave_wavelength_km"])
                  + float(fit["wave_phase_rad_at_pass_center"]))
    factor = float(fit["global_density_scale"]) * (
        1.0 + float(fit["linear_background_fraction_end_to_end_halfspan"])
        * x_norm[:, None, None]
        + float(fit["wave_amplitude_fraction_at_local_f2_peak"])
        * wave[:, None, None] * envelope)
    density *= factor
    if np.any(density <= 0) or not np.all(np.isfinite(density)):
        raise ValueError("Candidate density must be finite and positive")
    grid = build_pyiri_grid_from_axes(when, latitudes, longitudes, altitudes,
                                     f107=120.0, ap_daily=8.0,
                                     d_region_model="none")
    collision = grid.collision_freq
    if (grid.electron_temp_k is not None and grid.ion_temp_k is not None
            and grid.neutral_species_cm3 is not None):
        collision = effective_collision_frequency(
            grid.electron_temp_k, grid.ion_temp_k,
            density * 1e6, grid.neutral_species_cm3)
    candidate_grid = replace(grid, iono_en_grid=density,
                             iono_en_grid_5=density, collision_freq=collision)
    grid_output.parent.mkdir(parents=True, exist_ok=True)
    density_output.parent.mkdir(parents=True, exist_ok=True)
    save_ionosphere_grid_netcdf(grid_output, candidate_grid)
    np.savez_compressed(density_output, latitudes_deg=latitudes,
                        longitudes_deg=longitudes, altitudes_km=altitudes,
                        electron_density_cm3=density,
                        model=np.array("Ionogram-selected PyIRI latitude-wave fit"),
                        time_utc=np.array(when.isoformat()))
    print(json.dumps({"density": str(density_output), "forward_grid": str(grid_output),
                      "density_min_cm3": float(np.min(density)),
                      "density_max_cm3": float(np.max(density))}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit", type=Path, default=INITIAL_FIT)
    parser.add_argument("--density-output", type=Path, default=DENSITY_OUTPUT)
    parser.add_argument("--grid-output", type=Path, default=GRID_OUTPUT)
    args = parser.parse_args()
    build(args.fit, args.density_output, args.grid_output)
