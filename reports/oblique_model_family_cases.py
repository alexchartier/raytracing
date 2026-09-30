"""Prepare independent density truths and a common flexible oblique starting fit.

The ``truth`` action reads each external density grid and writes only simulated
observables and ray grids. ``start`` and ``candidate`` read the saved vertical
ionograms, 800 km in-situ samples, and PyIRI support grid, never the truth
density. Full oblique ray traces and candidate selection are separate jobs.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from ionogram_metrics import Ionogram  # noqa: E402
from python_raytrace.general_field_inverse import LocalDensity  # noqa: E402
from python_raytrace.general_topside_inverse import fit_spline_ionogram  # noqa: E402
from python_raytrace.geometry import llh_to_ecef  # noqa: E402
from python_raytrace.grid import (  # noqa: E402
    build_pyiri_grid_from_axes, load_ionosphere_grid_netcdf,
    save_ionosphere_grid_netcdf,
)
from python_raytrace.monotone_topside import monotone_topside_grid  # noqa: E402

RUN = ROOT / "reports/data/oblique_model_families"
SOURCES = RUN / "sources"
CASES = ("nequick", "chapman", "iri2016")
LAT_TX = -53.67859907285381
LON_TX = 7.7400542554165135
ALT_SC = 800.0
SEPARATION_KM = 600.0
SUPPORT_EPOCH = dt.datetime(2010, 1, 1, 12, tzinfo=dt.timezone.utc)
SUPPORT_F107 = 120.0
SUPPORT_AP = 8.0


def source_path(case: str) -> Path:
    return SOURCES / f"{case}_density.npz"


def vertical_path(case: str) -> Path:
    return SOURCES / f"{case}_vertical.npz"


def case_root(case: str) -> Path:
    return RUN / case


def save_private_grid(path: Path, grid) -> None:
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    save_ionosphere_grid_netcdf(path, grid)
    path.chmod(0o600)


def save_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    path.write_text(json.dumps(value, indent=2) + "\n")
    path.chmod(0o600)


def receiver_latitude() -> float:
    transmitter = np.asarray(llh_to_ecef(LAT_TX, LON_TX, ALT_SC * 1000.0))

    def distance(offset: float) -> float:
        receiver = np.asarray(llh_to_ecef(LAT_TX + offset, LON_TX, ALT_SC * 1000.0))
        return np.linalg.norm(receiver - transmitter) / 1000.0 - SEPARATION_KM

    return float(LAT_TX + brentq(distance, 0.1, 10.0))


def support() -> None:
    with np.load(source_path("iri2016"), allow_pickle=False) as source:
        lat = np.asarray(source["latitudes_deg"], dtype=float)
        lon = np.asarray(source["longitudes_deg"], dtype=float)
        alt = np.asarray(source["altitudes_km"], dtype=float)
    grid = build_pyiri_grid_from_axes(
        SUPPORT_EPOCH, lat, lon, alt, f107=SUPPORT_F107,
        ap_daily=SUPPORT_AP, d_region_model="none")
    save_private_grid(RUN / "support_grid.nc", grid)
    # The density prior deliberately comes from another latitude/longitude
    # and solar activity than the IRI-2016 validation profile.
    i = int(np.argmin(abs(lat - (-68.0))))
    j = int(np.argmin(abs(lon - (-12.0))))
    prior_profile = np.asarray(grid.iono_en_grid[i, j], dtype=float)
    prior_density = np.broadcast_to(prior_profile, grid.iono_en_grid.shape).copy()
    save_private_grid(RUN / "different_conditions_prior_grid.nc", replace(
        grid, iono_en_grid=prior_density, iono_en_grid_5=prior_density.copy(),
        metadata={**grid.metadata, "density_source": "PyIRI F10.7=120 at 68S, 12W; replicated prior"}))
    save_json(RUN / "support_provenance.json", {
        "support_epoch_utc": SUPPORT_EPOCH.isoformat(),
        "support_f107": SUPPORT_F107,
        "support_ap_daily": SUPPORT_AP,
        "density_prior_source_lat_lon_deg": [-68.0, -12.0],
        "truth_iri2016_f107": 72.6999969,
        "truth_iri2016_sounder_lat_lon_deg": [LAT_TX, LON_TX],
        "shared_support_fields": "PyIRI magnetic/collision fields shared by truth and candidates",
    })


def truth(case: str) -> None:
    if case not in CASES:
        raise ValueError(case)
    root = case_root(case)
    support_grid = load_ionosphere_grid_netcdf(RUN / "support_grid.nc")
    with np.load(source_path(case), allow_pickle=False) as source:
        lat = np.asarray(source["latitudes_deg"], dtype=float)
        lon = np.asarray(source["longitudes_deg"], dtype=float)
        alt = np.asarray(source["altitudes_km"], dtype=float)
        density = np.asarray(source["electron_density_cm3"], dtype=float)
        source_model = str(source["model"])
        source_time = str(source["time_utc"])
    if (not np.array_equal(lat, support_grid.latitudes_deg)
            or not np.array_equal(lon, support_grid.longitudes_deg)
            or not np.array_equal(alt, support_grid.altitudes_km)
            or density.shape != support_grid.iono_en_grid.shape
            or not np.all(np.isfinite(density)) or np.any(density <= 0)):
        raise ValueError("Truth axes or density do not match support grid")
    save_private_grid(root / "truth_grid.nc", replace(
        support_grid, iono_en_grid=density, iono_en_grid_5=density.copy(),
        metadata={**support_grid.metadata, "density_source": source_model}))
    rx_lat = receiver_latitude()
    distance_km = np.r_[np.arange(0.0, SEPARATION_KM, 1.0), SEPARATION_KM]
    locations = np.column_stack((LAT_TX + (rx_lat - LAT_TX) *
                                 distance_km / SEPARATION_KM,
                                 np.full(len(distance_km), LON_TX)))
    measured = RegularGridInterpolator((lat, lon, alt), density,
                                       bounds_error=True)(np.column_stack((
                                           locations, np.full(len(locations), ALT_SC))))
    save_json(root / "insitu_1km.json", {
        "provenance": "synthetic spacecraft density sampled from independent truth model",
        "altitude_km": ALT_SC, "spacing_km": 1.0,
        "distance_km": distance_km.tolist(),
        "lat_lon_deg": locations.tolist(),
        "density_cm3": measured.tolist(),
        "relative_uncertainty": 0.005,
    })
    save_json(root / "manifest.json", {
        "description": "Independent model-profile oblique test",
        "density_model": source_model,
        "source_time": source_time,
        "source_density": str(source_path(case).relative_to(ROOT)),
        "horizontal_structure": "IRI-2016 retains spatial variation; NeQuick and Chapman profiles are replicated horizontally",
        "profiles": [{"index": 1, "latitude_deg": LAT_TX,
                      "longitude_deg": LON_TX, "altitude_km": ALT_SC}],
        "satellite_separation_km": SEPARATION_KM,
        "frequency_sweep_mhz": [2.0, 10.0, 0.1],
        "homing_tolerance_m": 1000.0,
        "vertical_ionogram_observation": str(vertical_path(case).relative_to(ROOT)),
    })
    print(json.dumps({"case": case, "source_model": source_model,
                      "track_samples": len(measured),
                      "truth_fof2_mhz_at_transmitter": float(
                          .00898 * np.sqrt(np.max(RegularGridInterpolator(
                              (lat, lon), density)([[LAT_TX, LON_TX]])[0])))}))


def local_observations(case: str) -> LocalDensity:
    samples = json.loads((case_root(case) / "insitu_1km.json").read_text())
    return LocalDensity(np.asarray(samples["lat_lon_deg"]),
                        samples["altitude_km"],
                        np.asarray(samples["density_cm3"]),
                        samples["relative_uncertainty"])


def start(case: str) -> None:
    if case not in CASES:
        raise ValueError(case)
    observed = Ionogram.read(vertical_path(case))
    local = local_observations(case)
    fit = fit_spline_ionogram(observed.records, observed.frequencies, ALT_SC,
                             regularize_tail=True)
    layer = fit.layer
    support_grid = load_ionosphere_grid_netcdf(RUN / "support_grid.nc")
    unanchored = np.broadcast_to(layer.density_cm3(support_grid.altitudes_km),
                                 support_grid.iono_en_grid.shape).copy()
    base = replace(support_grid, iono_en_grid=unanchored,
                   iono_en_grid_5=unanchored.copy(),
                   metadata={**support_grid.metadata,
                             "density_source": "model-agnostic monotone spline fitted to O-mode vertical returns"})
    save_private_grid(case_root(case) / "unanchored_grid.nc", base)
    peak = float(layer.density_cm3(layer.hmf2_km))
    at_spacecraft = float(layer.density_cm3(ALT_SC))
    span = np.log(peak / at_spacecraft)
    shape = []
    for fraction in (.25, .60):
        z = layer.hmf2_km + fraction * (ALT_SC - layer.hmf2_km)
        shape.append(float(np.log(float(layer.density_cm3(z)) / at_spacecraft) / span))
    q1, q2 = shape
    if not 0.0 < q2 < q1 < 1.0:
        raise ValueError("Vertical spline did not yield a monotone shape")
    latitudes = np.array([LAT_TX, receiver_latitude()])
    initial = monotone_topside_grid(base, local, latitudes,
                                    np.full(2, q1), np.full(2, q2), np.zeros(2))
    save_private_grid(case_root(case) / "start/grid.nc", initial)
    save_json(case_root(case) / "start_parameters.json", {
        "vertical_observation": str(vertical_path(case).relative_to(ROOT)),
        "uses_full_truth_density": False,
        "uses_800km_track_observation": True,
        "prior_support": str((RUN / "support_grid.nc").relative_to(ROOT)),
        "density_start": "model-agnostic monotone spline from vertical O returns",
        "fof2_mhz": layer.fof2_mhz,
        "hmf2_km": layer.hmf2_km,
        "vertical_surrogate_range_mae_km": fit.range_mae_km[1],
        "q1_at_quarter": q1, "q2_at_three_fifths": q2,
        "latitude_controls_deg": latitudes.tolist(),
        "800km_observation_count": len(local.electron_density_cm3),
    })
    print(json.dumps({"case": case, "fof2_mhz": layer.fof2_mhz,
                      "hmf2_km": layer.hmf2_km, "q1": q1, "q2": q2}))


def candidate(case: str, name: str, dq1: float, dq2: float, dh: float,
              q1_gradient: float = 0.0, q2_gradient: float = 0.0) -> None:
    if case not in CASES or not name.replace("_", "").isalnum():
        raise ValueError("Invalid case or candidate name")
    params = json.loads((case_root(case) / "start_parameters.json").read_text())
    latitudes = np.asarray(params["latitude_controls_deg"])
    q1 = params["q1_at_quarter"] + dq1 + np.array([-q1_gradient, q1_gradient]) / 2
    q2 = params["q2_at_three_fifths"] + dq2 + np.array([-q2_gradient, q2_gradient]) / 2
    if np.any(q2 <= 0.02) or np.any(q1 >= .98) or np.any(q1 <= q2 + .08):
        raise ValueError("Proposed shape is outside monotone bounds")
    base = load_ionosphere_grid_netcdf(case_root(case) / "unanchored_grid.nc")
    updated = monotone_topside_grid(base, local_observations(case), latitudes,
                                    q1, q2, np.full(2, dh))
    dest = case_root(case) / name
    save_private_grid(dest / "grid.nc", updated)
    save_json(dest / "parameters.json", {
        "source": "oblique candidate from vertical spline and 800 km track samples",
        "uses_full_truth_density": False,
        "dq1": dq1, "dq2": dq2, "dh_km": dh,
        "q1_gradient": q1_gradient, "q2_gradient": q2_gradient,
        "q1_controls": q1.tolist(), "q2_controls": q2.tolist(),
    })
    print(json.dumps({"case": case, "candidate": name,
                      "q1": q1.tolist(), "q2": q2.tolist(), "dh_km": dh}))


def main() -> None:
    os.umask(0o077)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("support", "truth", "start", "candidate"))
    parser.add_argument("--case", choices=CASES)
    parser.add_argument("--name")
    parser.add_argument("--dq1", type=float, default=0.0)
    parser.add_argument("--dq2", type=float, default=0.0)
    parser.add_argument("--dh", type=float, default=0.0)
    parser.add_argument("--q1-gradient", type=float, default=0.0)
    parser.add_argument("--q2-gradient", type=float, default=0.0)
    args = parser.parse_args()
    if args.action == "support":
        support()
    elif args.case is None:
        parser.error("--case is required")
    elif args.action == "truth":
        truth(args.case)
    elif args.action == "start":
        start(args.case)
    elif args.name is None:
        parser.error("--name is required for a candidate")
    else:
        candidate(args.case, args.name, args.dq1, args.dq2, args.dh,
                  args.q1_gradient, args.q2_gradient)


if __name__ == "__main__":
    main()
