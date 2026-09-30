"""Oracle-only expressivity check for the model-independent field basis.

This diagnostic intentionally opens density truth to test whether the basis
can represent it. Its fitted coefficients must never be treated as a
retrieval or entered into ionogram candidate selection.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import lsq_linear

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from python_raytrace.general_field_inverse import (  # noqa: E402
    LocalDensity, TrackLogDensityBasis, candidate_grid,
)
from python_raytrace.general_topside_inverse import ChapmanF2  # noqa: E402
from python_raytrace.grid import load_ionosphere_grid_netcdf  # noqa: E402

DATA = ROOT / "reports/data"
OUTPUT = DATA / "general_field_iri_wave_structured/capacity.json"


def run() -> None:
    prior = load_ionosphere_grid_netcdf(
        DATA / "lat_wave_doppler_peak_round3/wide_grid.nc")
    rows = json.loads((DATA / "lat_wave_pass_manifest.json").read_text())["profiles"]
    locations = np.array([(row["latitude_deg"], row["longitude_deg"])
                          for row in rows])
    basis = TrackLogDensityBasis.from_locations(
        locations, spatial_centers=12,
        altitude_offsets_km=(-80., -30., 0., 45., 100., 180., 300., 500.),
        altitude_widths_km=(45., 40., 40., 50., 65., 90., 130., 170.))
    latitude, longitude = np.meshgrid(prior.latitudes_deg,
                                      prior.longitudes_deg, indexing="ij")
    along_km = basis.frame.project(latitude, longitude)
    phase = 2.0 * np.pi * along_km / 900.0
    independent_chapman = np.empty_like(prior.iono_en_grid)
    for (i, j), _ in np.ndenumerate(latitude):
        layer = ChapmanF2(
            fof2_mhz=5.30 + .42 * np.sin(phase[i, j]),
            hmf2_km=265.0 + 16.0 * np.cos(phase[i, j] + .5),
            bottomside_scale_km=52.0,
            topside_scale_km=85.0 + 22.0 * np.sin(phase[i, j] - .4),
            tail_curvature=.25 * np.cos(phase[i, j] + .3),
        )
        independent_chapman[i, j] = layer.density_cm3(prior.altitudes_km)
    with np.load(DATA / "lat_wave_pass_truth_density.npz", allow_pickle=False) as source:
        independent_iri = np.asarray(source["electron_density_cm3"], dtype=float)
    if independent_iri.shape != prior.iono_en_grid.shape:
        raise ValueError("IRI truth and prior grid shapes differ")

    def sample(field: np.ndarray) -> np.ndarray:
        return np.asarray(RegularGridInterpolator(
            (prior.latitudes_deg, prior.longitudes_deg), field,
            bounds_error=True)(locations), dtype=float)

    altitude = np.asarray(prior.altitudes_km)
    mask = (altitude >= 220) & (altitude <= 800)
    fit_altitude = altitude[mask]
    prior_profiles = sample(prior.iono_en_grid)
    ncoeff = int(np.prod(basis.coefficient_shape))
    design = np.column_stack([
        sample(basis.correction(
            prior, np.eye(ncoeff)[column].reshape(basis.coefficient_shape)))[:, mask].ravel()
        for column in range(ncoeff)
    ])
    # The peak gets extra weight, but the upper tail remains in the fit.
    weight = np.tile(1 + 8 * np.exp(-.5 * ((fit_altitude - 260) / 45) ** 2),
                     len(locations))
    cases = {}
    for label, truth_grid in (("independent_iri_wave", independent_iri),
                              ("analytic_chapman_wave", independent_chapman)):
        truth = sample(truth_grid)
        target = np.log(truth[:, mask] / prior_profiles[:, mask]).ravel()
        coefficient = lsq_linear(
            np.vstack((design * np.sqrt(weight[:, None]),
                       np.sqrt(.015) * np.eye(ncoeff))),
            np.r_[target * np.sqrt(weight), np.zeros(ncoeff)],
            bounds=(-.8, .8), max_iter=200).x.reshape(basis.coefficient_shape)
        print(label, "maximum absolute basis correction",
              float(np.max(abs(basis.correction(prior, coefficient)))), flush=True)
        measurement = LocalDensity(locations, 800.0,
                                   truth[:, np.argmin(abs(altitude - 800.0))])
        candidate = candidate_grid(prior, basis, coefficient,
                                   local_density=measurement,
                                   maximum_log_correction=2.0)
        retrieved = sample(candidate.iono_en_grid)
        true_peak = truth.max(axis=1)

        def metrics(profiles: np.ndarray) -> dict:
            return {
                "peak_density_mae_percent": float(np.mean(abs(
                    100 * (profiles.max(axis=1) / true_peak - 1)))),
                "hmf2_mae_km": float(np.mean(abs(
                    altitude[np.argmax(profiles, axis=1)]
                    - altitude[np.argmax(truth, axis=1)]))),
                "density_nrmse_220_to_800_percent": float(100 * np.mean(
                    np.sqrt(np.mean((profiles[:, mask] - truth[:, mask]) ** 2,
                                    axis=1)) / true_peak)),
            }

        cases[label] = {
            "prior": metrics(prior_profiles),
            "oracle_basis": metrics(retrieved),
            "maximum_log_density_correction": float(np.max(abs(
                basis.correction(prior, coefficient)))),
            "maximum_local_density_error_percent": float(np.max(abs(
                100 * (retrieved[:, np.argmin(abs(altitude - 800.0))]
                       / measurement.electron_density_cm3 - 1)))),
        }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    OUTPUT.write_text(json.dumps({
        "purpose": "oracle basis capacity only; truth density used to fit coefficients",
        "not_an_ionogram_retrieval": True,
        "spatial_centers": 12,
        "vertical_offsets_km": basis.altitude_offsets_km.tolist(),
        "cases": cases,
    }, indent=2) + "\n")
    OUTPUT.chmod(0o600)
    print(json.dumps(cases, indent=2))


if __name__ == "__main__":
    run()
