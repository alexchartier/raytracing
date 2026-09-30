"""Prepare and rank full-ray pilots for a model-independent sounder inverse.

Typical cycle: ``build`` writes a local-density-anchored baseline and paired
log-density perturbations; trace those grids with the existing complete O/X
generator; ``score`` ranks their ionograms; ``propose`` writes one combined
trust-region candidate. Only a traced candidate can be selected. This module
never opens a truth density grid.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import netCDF4
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from python_raytrace.general_field_inverse import (  # noqa: E402
    LocalDensity, TrackLogDensityBasis, candidate_grid, sampled_local_density,
    trust_region_step,
)
from python_raytrace.grid import (  # noqa: E402
    load_ionosphere_grid_netcdf, save_ionosphere_grid_netcdf,
)
from ionogram_metrics import Ionogram, score  # noqa: E402


def _private_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    path.write_text(text)
    path.chmod(0o600)


def _manifest(manifest: Path) -> tuple[list[dict], np.ndarray]:
    rows = json.loads(manifest.read_text())["profiles"]
    if len({int(row["index"]) for row in rows}) != len(rows):
        raise ValueError("Duplicate profile index")
    locations = np.array([(row["latitude_deg"], row["longitude_deg"])
                          for row in rows], dtype=float)
    return rows, locations


def _measurements(path: Path, locations: np.ndarray) -> LocalDensity:
    item = json.loads(path.read_text())
    if "density_cm3_by_profile" in item:
        density = np.asarray(item["density_cm3_by_profile"], dtype=float)
        altitude = float(item["altitude_km"])
    elif len(locations) == 1:
        density = np.array([float(item["electron_density_cm3"])])
        altitude = float(item["spacecraft_altitude_km"])
    else:
        raise ValueError("Expected one density per manifest profile")
    uncertainty = float(item.get("relative_uncertainty", 0.005))
    return LocalDensity(locations, altitude, density, uncertainty)


def _load_grid(path: Path):
    if sys.platform == "darwin":
        with netCDF4.Dataset(path) as dataset:
            cells = np.prod([len(dataset.dimensions[name]) for name in
                             ("lat", "lon", "alt")])
        if cells > 2_000_000:
            raise RuntimeError("Large field construction belongs on a Cartman compute node")
    return load_ionosphere_grid_netcdf(path)


def _save_candidate(path: Path, grid) -> None:
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    save_ionosphere_grid_netcdf(path, grid)
    path.chmod(0o600)


def _structured_directions(basis: TrackLogDensityBasis,
                           observed_dir: Path, rows: list[dict]) -> tuple[list[dict], dict]:
    """Four shared directions, including a wave inferred from O/X cutoffs."""
    spatial_count, vertical_count = basis.coefficient_shape
    if vertical_count < 5:
        raise ValueError("Structured pilots require five vertical basis functions")
    offsets = basis.altitude_offsets_km
    peak_index = int(np.argmin(abs(offsets)))
    low_index = int(np.argmin(abs(offsets + 50.0)))
    mid_index = int(np.argmin(abs(offsets - 80.0)))
    upper_index = int(np.argmin(abs(offsets - 210.0)))
    tail_index = int(np.argmin(abs(offsets - 420.0)))
    if len({peak_index, low_index, mid_index, upper_index, tail_index}) != 5:
        raise ValueError("Structured directions need separate low, peak, mid, upper, and tail modes")
    def direction(name: str, values: np.ndarray) -> dict:
        return {"name": name, "coefficients": values.tolist()}

    peak = np.zeros(basis.coefficient_shape)
    peak[:, peak_index] = 1.0
    height = np.zeros_like(peak)
    height[:, low_index], height[:, mid_index] = -1.0, 1.0
    width = np.zeros_like(peak)
    width[:, mid_index], width[:, upper_index] = -1.0, 1.0
    directions = [direction("peak", peak), direction("height", height),
                  direction("width", width)]
    sample_distance = basis.frame.sample_distance_km
    # Seven well-spaced samples still resolve wavelengths above twice the
    # largest along-track gap. The SAMI pass has seven observed positions.
    if len(rows) < 7 or np.ptp(sample_distance) < 300:
        tail = np.zeros_like(peak)
        tail[:, tail_index] = 1.0
        return directions + [direction("tail", tail)], {"wave_inferred": False}

    cutoff = []
    for row in rows:
        path = observed_dir / f"ionogram_{int(row['index']):02d}.npz"
        ionogram = Ionogram.read(path)
        # A Doppler gate selects near-vertical returns and can truncate the
        # observed nose severely on an irregular ray fan. Wave inference uses
        # the full accepted O/X cutoffs instead.
        valid = [ionogram.nose(mode) for mode in (1, -1)
                 if ionogram.nose(mode) is not None]
        cutoff.append(float(np.mean(valid)))
    cutoff = np.asarray(cutoff)
    s = sample_distance
    span = float(np.ptp(s))
    linear = np.column_stack((np.ones(len(s)), (s - np.mean(s)) / 1000.0))
    linear_fit = np.linalg.lstsq(linear, cutoff, rcond=None)[0]
    baseline_sse = float(np.mean((linear @ linear_fit - cutoff) ** 2))
    low = max(300.0, 2.0 * float(np.max(np.diff(np.sort(s)))), span / 7.0)
    high = min(2200.0, span * 1.5)
    candidates = []
    for wavelength in np.linspace(low, high, 220):
        phase = 2 * np.pi * s / wavelength
        design = np.column_stack((linear, np.cos(phase), np.sin(phase)))
        coefficients = np.linalg.lstsq(design, cutoff, rcond=None)[0]
        error = float(np.mean((design @ coefficients - cutoff) ** 2))
        candidates.append((error, wavelength, coefficients))
    error, wavelength, coefficients = min(candidates, key=lambda item: item[0])
    gain = ((baseline_sse - error) / baseline_sse if baseline_sse > 0 else 0.0)
    if gain < .12:
        spatial_pattern = ((basis.center_distance_km - np.mean(s))
                           / max(span / 2.0, 1.0))
        label = "gradient"
    else:
        phase_centers = 2 * np.pi * basis.center_distance_km / wavelength
        spatial_pattern = (coefficients[2] * np.cos(phase_centers)
                           + coefficients[3] * np.sin(phase_centers))
        label = "wave"
    spatial_pattern /= max(np.max(abs(spatial_pattern)), 1e-12)
    structured = np.zeros_like(peak)
    structured[:, peak_index] = spatial_pattern
    directions.append(direction(label, structured))
    return directions, {"wave_inferred": label == "wave",
                        "wavelength_km": float(wavelength),
                        "harmonic_amplitude_mhz": float(np.linalg.norm(coefficients[2:])),
                        "variance_reduction_over_linear": float(gain),
                        "cutoff_mhz_by_profile": cutoff.tolist()}


def build(args: argparse.Namespace) -> None:
    os.umask(0o077)
    rows, locations = _manifest(args.manifest)
    if (args.perturbation <= 0 or args.maximum_log_correction <= 0
            or any(index not in {int(row["index"]) for row in rows}
                   for index in args.pilot_indices)):
        raise ValueError("Invalid perturbation, correction bound, or pilot profile")
    measurement = _measurements(args.local_density, locations)
    basis_options = {}
    if args.altitude_offsets_km is not None or args.altitude_widths_km is not None:
        if args.altitude_offsets_km is None or args.altitude_widths_km is None:
            raise ValueError("Supply both altitude offsets and widths")
        basis_options = {"altitude_offsets_km": tuple(args.altitude_offsets_km),
                         "altitude_widths_km": tuple(args.altitude_widths_km)}
    basis = TrackLogDensityBasis.from_locations(
        locations, spatial_centers=args.spatial_centers, **basis_options)
    prior = _load_grid(args.prior_grid)
    modes = sorted(set(args.vertical_modes))
    if not modes or min(modes) < 0 or max(modes) >= basis.coefficient_shape[1]:
        raise ValueError("Vertical modes must index the configured basis")
    active_spatial = (list(range(basis.coefficient_shape[0]))
                      if args.active_spatial_centers is None else
                      sorted(set(args.active_spatial_centers)))
    if not active_spatial or min(active_spatial) < 0 or max(active_spatial) >= basis.coefficient_shape[0]:
        raise ValueError("Active spatial centers must index the configured basis")
    zero = np.zeros(basis.coefficient_shape)
    candidates = []
    directions = []
    wave_diagnostic = None

    def write(name: str, c: np.ndarray) -> None:
        trial = candidate_grid(prior, basis, c, local_density=measurement,
                               maximum_log_correction=args.maximum_log_correction)
        output = args.output / name / "grid.nc"
        _save_candidate(output, trial)
        local_error = sampled_local_density(trial, measurement)
        candidates.append({"name": name, "grid": str(output),
                           "coefficients": c.tolist(),
                           "maximum_local_density_error_percent": float(np.max(abs(
                               100 * (local_error / measurement.electron_density_cm3 - 1))))})

    write("baseline", zero)
    if args.pilot_design == "structured":
        directions, wave_diagnostic = _structured_directions(
            basis, args.observed_dir, rows)
        for row in directions:
            coefficients = np.asarray(row["coefficients"], dtype=float)
            if args.active_spatial_centers is not None:
                inactive = sorted(set(range(basis.coefficient_shape[0]))
                                  - set(active_spatial))
                coefficients[inactive] = 0.0
                row["coefficients"] = coefficients.tolist()
            for sign, prefix in ((1, "plus"), (-1, "minus")):
                write(f"{prefix}_{row['name']}", sign * args.perturbation * coefficients)
    else:
        for spatial in active_spatial:
            for vertical in modes:
                for sign, prefix in ((1, "plus"), (-1, "minus")):
                    coefficients = zero.copy()
                    coefficients[spatial, vertical] = sign * args.perturbation
                    write(f"{prefix}_{spatial:02d}_{vertical:02d}", coefficients)
    result = {
        "selection_uses_truth_density": False,
        "prior_grid": str(args.prior_grid),
        "manifest": str(args.manifest),
        "observed_ionograms": str(args.observed_dir),
        "local_density_observations": str(args.local_density),
        "profile_indices": [int(row["index"]) for row in rows],
        "pilot_profile_indices": args.pilot_indices,
        "spatial_centers": args.spatial_centers,
        "active_spatial_centers": active_spatial,
        "active_vertical_modes": modes,
        "spatial_center_distance_km": basis.center_distance_km.tolist(),
        "spatial_width_km": basis.spatial_width_km,
        "vertical_offsets_km": basis.altitude_offsets_km.tolist(),
        "vertical_widths_km": basis.altitude_widths_km.tolist(),
        "perturbation": args.perturbation,
        "pilot_design": args.pilot_design,
        "directions": directions,
        "wave_diagnostic": wave_diagnostic,
        "maximum_log_correction": args.maximum_log_correction,
        "candidate_ionogram_requirement":
            "complete 2-10 MHz O/X sweep, same 100 kHz grid and gap-recovery settings as observations",
        "candidates": candidates,
    }
    _private_write(args.output / "plan.json", json.dumps(result, indent=2) + "\n")
    print(json.dumps({"candidate_count": len(candidates),
                      "pilot_profile_indices": args.pilot_indices,
                      "maximum_local_density_error_percent": max(
                          row["maximum_local_density_error_percent"]
                          for row in candidates)}, indent=2))


def _gated_nose(path: Path, mode: int, gate_hz: float) -> float | None:
    with np.load(path, allow_pickle=False) as source:
        if "spacecraft_doppler_hz" not in source:
            return None
        records = np.asarray(source["records"], dtype=float)
        doppler = np.rint(np.asarray(source["spacecraft_doppler_hz"], dtype=float))
        frequencies = np.asarray(source["frequencies_mhz"], dtype=float)
    if len(doppler) != len(records):
        raise ValueError(f"Doppler and return arrays differ in {path}")
    selected = records[(records[:, 1] == mode) & (abs(doppler) <= gate_hz)]
    return float(frequencies[int(np.max(selected[:, 0]))]) if len(selected) else None


def score_pair(observed_path: Path, predicted_path: Path,
               gate_hz: float = 15.0) -> dict:
    observed = Ionogram.read(observed_path)
    predicted = Ionogram.read(predicted_path)
    standard = score(observed, predicted)
    with np.load(observed_path, allow_pickle=False) as source:
        observed_doppler = "spacecraft_doppler_hz" in source
    with np.load(predicted_path, allow_pickle=False) as source:
        predicted_doppler = "spacecraft_doppler_hz" in source
    if gate_hz <= 0 or not (observed_doppler and predicted_doppler):
        return {**standard, "combined": standard["total"], "gated_cutoff": None}
    gate_errors = []
    for mode in (1, -1):
        a = _gated_nose(observed_path, mode, gate_hz)
        b = _gated_nose(predicted_path, mode, gate_hz)
        gate_errors.append(1.0 if a is None or b is None
                           else min(abs(a - b) / .8, 1.0))
    gated = float(np.mean(gate_errors))
    return {**standard,
            "combined": .7 * standard["total"] + .3 * gated,
            "gated_cutoff": gated}


def score_pilots(args: argparse.Namespace) -> None:
    plan = json.loads(args.plan.read_text())
    indices = (args.indices if args.indices is not None else
               plan["pilot_profile_indices"])
    if not indices:
        raise ValueError("No pilot profiles")
    rows = []
    observed = Path(plan["observed_ionograms"])
    known_names = {row["name"] for row in plan["candidates"]}
    candidate_names = ([row["name"] for row in plan["candidates"]]
                       if args.candidate_names is None else args.candidate_names)
    if any(name not in known_names for name in candidate_names):
        raise ValueError("Requested candidate is not in the plan")
    candidate_names += args.extra_candidate_names
    reference_dirs = {}
    for item in args.reference:
        name, separator, directory = item.partition("=")
        if not separator or not name or not directory:
            raise ValueError("Reference must be NAME=IONOGRAM_DIRECTORY")
        reference_dirs[name] = Path(directory)
    candidate_names += list(reference_dirs)
    if len(candidate_names) != len(set(candidate_names)):
        raise ValueError("Duplicate candidate name")
    for name in candidate_names:
        forward = reference_dirs.get(name, args.forward_root / name / "ionograms")
        parts = []
        for index in indices:
            filename = f"ionogram_{index:02d}.npz"
            parts.append(score_pair(observed / filename, forward / filename,
                                    args.doppler_gate_hz))
        rows.append({"name": name, "mean_combined_score": float(np.mean(
                         [part["combined"] for part in parts])),
                     "mean_standard_score": float(np.mean(
                         [part["total"] for part in parts])),
                     "by_profile": parts})
    result = {"selection_uses_truth_density": False,
              "plan": str(args.plan), "profile_indices": indices,
              "scores": rows}
    _private_write(args.output, json.dumps(result, indent=2) + "\n")
    print(json.dumps({"best": min(rows, key=lambda row: row["mean_combined_score"])["name"],
                      "scores": {row["name"]: row["mean_combined_score"]
                                 for row in rows}}, indent=2))


def propose(args: argparse.Namespace) -> None:
    plan = json.loads(args.plan.read_text())
    scored = json.loads(args.scores.read_text())
    if scored["plan"] != str(args.plan):
        raise ValueError("Pilot scores refer to a different plan")
    by_name = {row["name"]: row["mean_combined_score"] for row in scored["scores"]}
    if plan["pilot_design"] == "structured":
        labels = [row["name"] for row in plan["directions"]]
        plus = np.array([by_name[f"plus_{name}"] for name in labels])
        minus = np.array([by_name[f"minus_{name}"] for name in labels])
    else:
        active = [(spatial, vertical)
                  for spatial in plan["active_spatial_centers"]
                  for vertical in plan["active_vertical_modes"]]
        plus = np.array([by_name[f"plus_{i:02d}_{j:02d}"] for i, j in active])
        minus = np.array([by_name[f"minus_{i:02d}_{j:02d}"] for i, j in active])
    step = trust_region_step(by_name["baseline"], plus, minus,
                             perturbation=plan["perturbation"],
                             radius=args.radius, regularization=args.regularization)
    coefficients = np.zeros((plan["spatial_centers"],
                             len(plan["vertical_offsets_km"])))
    if plan["pilot_design"] == "structured":
        for row, value in zip(plan["directions"], step):
            coefficients += value * np.asarray(row["coefficients"], dtype=float)
    else:
        for (spatial, vertical), value in zip(active, step):
            coefficients[spatial, vertical] = value
    _, locations = _manifest(Path(plan["manifest"]))
    basis = TrackLogDensityBasis.from_locations(
        locations, spatial_centers=plan["spatial_centers"],
        altitude_offsets_km=tuple(plan["vertical_offsets_km"]),
        altitude_widths_km=tuple(plan["vertical_widths_km"]))
    measurement = _measurements(Path(plan["local_density_observations"]), locations)
    prior = _load_grid(Path(plan["prior_grid"]))
    candidate = candidate_grid(prior, basis, coefficients, local_density=measurement,
                               maximum_log_correction=plan["maximum_log_correction"])
    output = args.output / "grid.nc"
    _save_candidate(output, candidate)
    _private_write(args.output / "proposal.json", json.dumps({
        "selection_uses_truth_density": False,
        "full_ray_ionogram_selection_complete": False,
        "source_plan": str(args.plan), "pilot_scores": str(args.scores),
        "coefficients": coefficients.tolist(),
        "grid": str(output),
        "pilot_best_single_candidate": min(scored["scores"],
             key=lambda row: row["mean_combined_score"])["name"],
        "decision": "trace combined grid on all profiles, then compare with baseline and best pilot",
    }, indent=2) + "\n")
    print(str(output))


def select_full(args: argparse.Namespace) -> None:
    """Freeze an all-profile choice before opening any density truth."""
    scored = json.loads(args.scores.read_text())
    plan = json.loads(Path(scored["plan"]).read_text())
    if sorted(scored["profile_indices"]) != sorted(plan["profile_indices"]):
        raise ValueError("Final selection requires every profile in the manifest")
    by_name = {row["name"]: row for row in scored["scores"]}
    if args.reference_name not in by_name:
        raise ValueError("The reference ionograms must be among the scored candidates")
    reference = by_name[args.reference_name]
    winner = min(scored["scores"], key=lambda row: row["mean_combined_score"])
    differences = np.array([a["combined"] - b["combined"]
                            for a, b in zip(winner["by_profile"],
                                            reference["by_profile"])])
    rng = np.random.default_rng(20260930)
    bootstrap = rng.choice(differences, size=(10_000, len(differences)),
                           replace=True).mean(axis=1)
    result = {
        "selection_uses_truth_density": False,
        "chosen_candidate": winner["name"],
        "reference_candidate": args.reference_name,
        "profile_indices": scored["profile_indices"],
        "chosen_mean_combined_score": winner["mean_combined_score"],
        "reference_mean_combined_score": reference["mean_combined_score"],
        "paired_score_change": float(np.mean(differences)),
        "positions_with_better_score": int(np.count_nonzero(differences < 0)),
        "paired_position_bootstrap_95_percent_interval": np.quantile(
            bootstrap, [.025, .975]).tolist(),
        "full_ray_scores": str(args.scores),
    }
    _private_write(args.output, json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    p = commands.add_parser("build")
    p.add_argument("--prior-grid", type=Path, required=True)
    p.add_argument("--manifest", type=Path, required=True)
    p.add_argument("--observed-dir", type=Path, required=True)
    p.add_argument("--local-density", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--spatial-centers", type=int, default=4)
    p.add_argument("--altitude-offsets-km", nargs="+", type=float)
    p.add_argument("--altitude-widths-km", nargs="+", type=float)
    p.add_argument("--pilot-design", choices=("structured", "coordinate"),
                   default="structured")
    p.add_argument("--vertical-modes", nargs="+", type=int, default=[1, 2, 3, 4])
    p.add_argument("--active-spatial-centers", nargs="+", type=int)
    p.add_argument("--pilot-indices", nargs="+", type=int, default=[3, 8, 14, 18])
    p.add_argument("--perturbation", type=float, default=.05)
    p.add_argument("--maximum-log-correction", type=float, default=.5)
    p.set_defaults(func=build)
    p = commands.add_parser("score")
    p.add_argument("--plan", type=Path, required=True)
    p.add_argument("--forward-root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--indices", nargs="+", type=int)
    p.add_argument("--candidate-names", nargs="*")
    p.add_argument("--doppler-gate-hz", type=float, default=15.0)
    p.add_argument("--extra-candidate-names", nargs="+", default=[])
    p.add_argument("--reference", action="append", default=[],
                   help="NAME=directory for a frozen full-ray candidate")
    p.set_defaults(func=score_pilots)
    p = commands.add_parser("propose")
    p.add_argument("--plan", type=Path, required=True)
    p.add_argument("--scores", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--radius", type=float, default=.12)
    p.add_argument("--regularization", type=float, default=1.0)
    p.set_defaults(func=propose)
    p = commands.add_parser("select")
    p.add_argument("--scores", type=Path, required=True)
    p.add_argument("--reference-name", default="previous")
    p.add_argument("--output", type=Path, required=True)
    p.set_defaults(func=select_full)
    parsed = parser.parse_args()
    parsed.func(parsed)
