"""Ionogram-only pilot inversion of the SAMI3 wave, followed by held-out evaluation."""

from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from python_raytrace.grid import load_ionosphere_grid_netcdf, save_ionosphere_grid_netcdf
from ionogram_metrics import Ionogram, score

RUN = ROOT / "reports/data/sami3_wave_20170111_0600"
FIGURES = ROOT / "reports/figures"
INDICES = (1, 5, 9, 12, 14, 17, 20)


def _profiles() -> list[dict]:
    return json.loads((RUN / "manifest.json").read_text())["profiles"]


def _ridge_residual(obs: Ionogram, prior: Ionogram) -> float:
    """Median low-frequency O-mode group-range residual in km."""
    a, b = obs.ridge(1), prior.ridge(1)
    common = sorted(set(a) & set(b))
    values = [a[i] - b[i] for i in common if 2.5 <= obs.frequencies[i] <= 3.7]
    return float(np.median(values)) if len(values) >= 3 else float("nan")


def _crossing_height(alt: np.ndarray, profile: np.ndarray, frequency: float) -> float:
    peak = int(np.argmax(profile))
    topside = profile[peak:]
    target = (frequency / .00898) ** 2
    if target >= topside[0] or target <= topside[-1]:
        return float("nan")
    return float(np.interp(target, topside[::-1], alt[peak:][::-1]))


def _peak_density_and_height(alt: np.ndarray, profiles: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Three-point quadratic peak on the shared, uniform altitude grid."""
    rows = np.arange(len(profiles))
    index = np.clip(np.argmax(profiles, axis=1), 1, len(alt) - 2)
    left = profiles[rows, index - 1]
    center = profiles[rows, index]
    right = profiles[rows, index + 1]
    curvature = left - 2 * center + right
    offset = np.zeros_like(center)
    np.divide(.5 * (left - right), curvature, out=offset, where=curvature < 0)
    offset = np.clip(offset, -1.0, 1.0)
    height = alt[index] + offset * (alt[1] - alt[0])
    density = center - .25 * (left - right) * offset
    return density, height


def infer_observables() -> dict:
    grid = load_ionosphere_grid_netcdf(RUN / "pyiri_prior_grid.nc")
    rows = []
    for index in INDICES:
        row = _profiles()[index - 1]
        observed = Ionogram.read(RUN / "vertical_truth" / f"ionogram_{index:02d}.npz")
        prior = Ionogram.read(RUN / "vertical_prior" / f"ionogram_{index:02d}.npz")
        lat = float(row["latitude_deg"])
        o_nose, p_nose = observed.nose(1), prior.nose(1)
        x_nose, px_nose = observed.nose(-1), prior.nose(-1)
        pairs = [(o_nose, p_nose), (x_nose, px_nose)]
        ratios = [(o / p) ** 2 for o, p in pairs if o and p and o < 9.9 and p < 9.9]
        if not ratios:
            raise ValueError(f"No useful O/X nose at profile {index}")
        peak_factor = float(np.exp(np.mean(np.log(ratios))))
        p = RegularGridInterpolator((grid.latitudes_deg, grid.longitudes_deg),
                                    grid.iono_en_grid)([[lat, float(row["longitude_deg"])]])[0]
        h0 = float(grid.altitudes_km[np.argmax(p)])
        low_residual = _ridge_residual(observed, prior)
        h_cross = _crossing_height(grid.altitudes_km, p, 3.2)
        # A broader topside generally raises dispersive group range. The
        # 250 km per unit stretch is only a proposal; full rays select gain.
        stretch = 1.0 + low_residual / 250.0
        stretch = float(np.clip(stretch, .65, 2.7)) if np.isfinite(stretch) else 1.0
        rows.append({"index": index, "latitude_deg": lat,
                     "o_nose_mhz": o_nose, "x_nose_mhz": x_nose,
                     "prior_o_nose_mhz": p_nose, "prior_x_nose_mhz": px_nose,
                     "peak_density_factor": peak_factor,
                     "low_frequency_range_residual_km": low_residual,
                     "prior_3p2mhz_crossing_km": h_cross,
                     "proposed_topside_stretch": stretch,
                     "prior_score": score(observed, prior)})
    result = {"indices": list(INDICES), "selection_uses_truth_density": False,
              "method": "O/X nose ratio for peak scaling; low-frequency O group-range residual for topside stretch proposal; all candidates forward checked",
              "profiles": rows}
    (RUN / "vertical_fit.json").write_text(json.dumps({**result, "wave_wavelength_km": 1300.0}, indent=2) + "\n")
    return result


def build_candidate(name: str, stretch_gain: float) -> dict:
    fit = json.loads((RUN / "vertical_fit.json").read_text())
    grid = load_ionosphere_grid_netcdf(RUN / "pyiri_prior_grid.nc")
    sampled_lat = np.array([r["latitude_deg"] for r in fit["profiles"]])
    peak_factors = np.array([r["peak_density_factor"] for r in fit["profiles"]])
    proposed_stretches = np.array([r["proposed_topside_stretch"] for r in fit["profiles"]])
    lat = np.asarray(grid.latitudes_deg)
    peak_factor = np.clip(PchipInterpolator(sampled_lat, peak_factors,
                                            extrapolate=False)(np.clip(lat, sampled_lat[0], sampled_lat[-1])), .65, 1.5)
    stretch = np.clip(1.0 + stretch_gain * (PchipInterpolator(
        sampled_lat, proposed_stretches, extrapolate=False)(
            np.clip(lat, sampled_lat[0], sampled_lat[-1])) - 1.0), .7, 2.7)
    alt = np.asarray(grid.altitudes_km)
    base = np.asarray(grid.iono_en_grid)
    revised = np.empty_like(base)
    for i in range(len(lat)):
        for j in range(len(grid.longitudes_deg)):
            profile = base[i, j]
            h = alt[np.argmax(profile)]
            sampled = np.where(alt > h, h + (alt - h) / stretch[i], alt)
            revised[i, j] = np.interp(sampled, alt, profile) * peak_factor[i]
    directory = RUN / "vertical_candidates" / name
    directory.mkdir(parents=True, exist_ok=True)
    save_ionosphere_grid_netcdf(directory / "grid.nc",
                                replace(grid, iono_en_grid=revised, iono_en_grid_5=revised.copy(),
                                        metadata={**grid.metadata, "retrieval": name}))
    np.savez_compressed(directory / "density.npz", latitudes_deg=grid.latitudes_deg,
                        longitudes_deg=grid.longitudes_deg,
                        altitudes_km=grid.altitudes_km, electron_density_cm3=revised,
                        model="PyIRI plus ionogram-inferred peak/topside correction")
    summary = {"name": name, "stretch_gain": stretch_gain,
               "peak_factor_range": [float(peak_factor.min()), float(peak_factor.max())],
               "stretch_range": [float(stretch.min()), float(stretch.max())]}
    (directory / "candidate.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def refine_nose(name: str, start_name: str, gain: float = 1.0,
                sigma_km: float = 60.0) -> dict:
    """Use remaining O/X nose errors to update the peak without truth density."""
    if not 0 < gain <= 1.5 or not 20 <= sigma_km <= 150:
        raise ValueError("Invalid nose-refinement gain or width")
    start = RUN / "vertical_candidates" / start_name
    grid = load_ionosphere_grid_netcdf(start / "grid.nc")
    rows = []
    for index in INDICES:
        observed = Ionogram.read(RUN / "vertical_truth" / f"ionogram_{index:02d}.npz")
        modeled = Ionogram.read(start / "ionograms" / f"ionogram_{index:02d}.npz")
        pairs = [(observed.nose(mode), modeled.nose(mode)) for mode in (1, -1)]
        ratios = [(a / b) ** 2 for a, b in pairs if a and b and a < 9.9 and b < 9.9]
        if not ratios:
            raise ValueError(f"No paired O/X nose at profile {index}")
        rows.append({"index": index, "latitude_deg": float(_profiles()[index - 1]["latitude_deg"]),
                     "observed_noses_mhz": [a for a, _ in pairs],
                     "modeled_noses_mhz": [b for _, b in pairs],
                     "peak_factor_proposal": float(np.exp(np.mean(np.log(ratios))))})
    sampled_lat = np.array([r["latitude_deg"] for r in rows])
    factors = np.array([r["peak_factor_proposal"] for r in rows])
    lat = np.asarray(grid.latitudes_deg)
    applied = np.clip(np.exp(gain * PchipInterpolator(sampled_lat, np.log(factors))(
        np.clip(lat, sampled_lat[0], sampled_lat[-1]))), .8, 1.2)
    density = np.asarray(grid.iono_en_grid)
    alt = np.asarray(grid.altitudes_km)
    h = alt[np.argmax(density, axis=2)]
    envelope = np.exp(-.5 * ((alt[None, None, :] - h[:, :, None]) / sigma_km) ** 2)
    revised = density * (1 + (applied[:, None, None] - 1) * envelope)
    directory = RUN / "vertical_candidates" / name
    directory.mkdir(parents=True, exist_ok=True)
    save_ionosphere_grid_netcdf(directory / "grid.nc",
                                replace(grid, iono_en_grid=revised, iono_en_grid_5=revised.copy(),
                                        metadata={**grid.metadata, "retrieval": name}))
    np.savez_compressed(directory / "density.npz", latitudes_deg=grid.latitudes_deg,
                        longitudes_deg=grid.longitudes_deg,
                        altitudes_km=grid.altitudes_km, electron_density_cm3=revised,
                        model="PyIRI ionogram-inferred second-round peak correction")
    result = {"name": name, "start_name": start_name, "gain": gain,
              "sigma_km": sigma_km, "selection_uses_truth_density": False,
              "applied_peak_factor_range": [float(applied.min()), float(applied.max())],
              "profiles": rows}
    (directory / "candidate.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def select(names: list[str]) -> dict:
    rows = []
    for name in ["prior", *names]:
        modeled_dir = (RUN / "vertical_prior" if name == "prior" else
                       RUN / "vertical_candidates" / name / "ionograms")
        per_profile = []
        for index in INDICES:
            obs = Ionogram.read(RUN / "vertical_truth" / f"ionogram_{index:02d}.npz")
            modeled = Ionogram.read(modeled_dir / f"ionogram_{index:02d}.npz")
            per_profile.append({"index": index, **score(obs, modeled)})
        rows.append({"name": name,
                     "mean_total_score": float(np.mean([v["total"] for v in per_profile])),
                     "profiles": per_profile})
    selected = min(rows, key=lambda r: r["mean_total_score"])["name"]
    result = {"selected": selected, "selection_uses_truth_density": False,
              "indices": list(INDICES), "candidates": rows}
    (RUN / "vertical_selection.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def evaluate() -> dict:
    selection = json.loads((RUN / "vertical_selection.json").read_text())
    name = selection["selected"]
    candidate = (RUN / "pyiri_prior_grid.nc" if name == "prior" else
                 RUN / "vertical_candidates" / name / "grid.nc")
    truth = np.asarray(load_ionosphere_grid_netcdf(RUN / "truth_grid.nc").iono_en_grid)
    model_grid = load_ionosphere_grid_netcdf(candidate)
    modeled = np.asarray(model_grid.iono_en_grid)
    lat, lon, alt = model_grid.latitudes_deg, model_grid.longitudes_deg, model_grid.altitudes_km
    positions = np.array([_profiles()[i - 1]["latitude_deg"] for i in INDICES])
    points = np.column_stack((positions, np.full(len(positions), 190.0)))
    a = RegularGridInterpolator((lat, lon), truth)(points)
    b = RegularGridInterpolator((lat, lon), modeled)(points)
    pa, ha = _peak_density_and_height(alt, a)
    pb, hb = _peak_density_and_height(alt, b)
    mask = (alt >= 150) & (alt <= 600)
    result = {"selected": name, "positions_deg": positions.tolist(),
              "truth_fof2_mhz": (.00898 * np.sqrt(pa)).tolist(),
              "retrieved_fof2_mhz": (.00898 * np.sqrt(pb)).tolist(),
              "truth_hmf2_km": ha.tolist(), "retrieved_hmf2_km": hb.tolist(),
              "mean_absolute_fof2_error_mhz": float(np.mean(abs(.00898 * np.sqrt(pb) - .00898 * np.sqrt(pa)))),
              "mean_absolute_peak_density_error_percent": float(np.mean(abs(pb / pa - 1)) * 100),
              "mean_absolute_hmf2_error_km": float(np.mean(abs(hb - ha))),
              "normalized_density_rms_percent": float(100 * np.sqrt(np.mean((a[:, mask] - b[:, mask]) ** 2)) / np.max(a[:, mask]))}
    (RUN / "vertical_evaluation.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def evaluate_oblique(name: str) -> dict:
    """Post-selection density and ionogram evaluation at the seven link midpoints."""
    directory = RUN / "oblique_wave_600km"
    os.environ["SOUNDER_CASE_ROOT"] = str(RUN)
    from oblique_wave_pass import doppler_score
    modeled_ions = directory / ("baseline" if name == "baseline" else f"{name}/ionograms")
    model_path = (RUN / "vertical_selected_grid.nc" if name == "baseline" else
                  directory / name / "grid.nc")
    truth_grid = load_ionosphere_grid_netcdf(RUN / "truth_grid.nc")
    model_grid = load_ionosphere_grid_netcdf(model_path)
    midpoints = []
    scores = []
    doppler_scores = []
    for index in INDICES:
        source = directory / "truth" / f"ionogram_{index:02d}.npz"
        target = modeled_ions / f"ionogram_{index:02d}.npz"
        with np.load(source, allow_pickle=False) as d:
            midpoints.append([.5 * (float(d["tx_lat_deg"]) + float(d["rx_lat_deg"])),
                              float(d["tx_lon_deg"])])
        scores.append(score(Ionogram.read(source), Ionogram.read(target)))
        doppler_scores.append(doppler_score(source, target))
    points = np.array(midpoints)
    axes = (truth_grid.latitudes_deg, truth_grid.longitudes_deg)
    truth = RegularGridInterpolator(axes, truth_grid.iono_en_grid)(points)
    modeled = RegularGridInterpolator(axes, model_grid.iono_en_grid)(points)
    altitude = truth_grid.altitudes_km
    a, h_a = _peak_density_and_height(altitude, truth)
    b, h_b = _peak_density_and_height(altitude, modeled)
    mask = (altitude >= 150) & (altitude <= 600)
    result = {"selected": name, "link_midpoint_latitudes_deg": points[:, 0].tolist(),
              "mean_ionogram_score": float(np.mean([row["total"] for row in scores])),
              "mean_1hz_doppler_score": float(np.mean(doppler_scores)),
              "ionogram_scores": scores,
              "doppler_scores": doppler_scores,
              "truth_fof2_mhz": (.00898 * np.sqrt(a)).tolist(),
              "retrieved_fof2_mhz": (.00898 * np.sqrt(b)).tolist(),
              "truth_hmf2_km": h_a.tolist(), "retrieved_hmf2_km": h_b.tolist(),
              "mean_absolute_fof2_error_mhz": float(np.mean(abs(.00898 * np.sqrt(b) - .00898 * np.sqrt(a)))),
              "mean_absolute_peak_density_error_percent": float(np.mean(abs(b / a - 1)) * 100),
              "mean_absolute_hmf2_error_km": float(np.mean(abs(h_b - h_a))),
              "normalized_density_rms_percent": float(100 * np.sqrt(np.mean((truth[:, mask] - modeled[:, mask]) ** 2)) / np.max(truth[:, mask]))}
    (RUN / "oblique_evaluation.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("infer")
    b = sub.add_parser("build")
    b.add_argument("--name", required=True)
    b.add_argument("--stretch-gain", type=float, required=True)
    n = sub.add_parser("refine-nose")
    n.add_argument("--name", required=True)
    n.add_argument("--start", required=True)
    n.add_argument("--gain", type=float, default=1.0)
    n.add_argument("--sigma-km", type=float, default=60.0)
    s = sub.add_parser("select")
    s.add_argument("--names", nargs="+", required=True)
    sub.add_parser("evaluate")
    o = sub.add_parser("evaluate-oblique")
    o.add_argument("--name", required=True)
    args = parser.parse_args()
    if args.command == "infer": result = infer_observables()
    elif args.command == "build": result = build_candidate(args.name, args.stretch_gain)
    elif args.command == "refine-nose": result = refine_nose(args.name, args.start,
                                                               args.gain, args.sigma_km)
    elif args.command == "select": result = select(args.names)
    elif args.command == "evaluate": result = evaluate()
    else: result = evaluate_oblique(args.name)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
