"""Fit a NeQuick-G ionogram using only its saved O/X returns.

The retrieval never imports NeQuick. For the local-density case, preparation
exposes one scalar at the spacecraft. The remaining independent density and
its peak parameters are opened only by ``evaluate``, after full-ray selection.
The default July profile became a development case when its O/X cutoff
difference was inspected; do not label it a prospective held-out validation
of that rule. A separate case may be supplied with ``--run``.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "reports"))

from ionogram_metrics import Ionogram, score  # noqa: E402
from recover_missing_ionogram_returns import interior_gaps  # noqa: E402
from python_raytrace.general_topside_inverse import (  # noqa: E402
    PLASMA_MHZ_PER_SQRT_CM3, MonotoneSplineF2, fit_spline_ionogram,
)

RUN = ROOT / "reports/data/nequick_heldout_2026_case"
OLD = "legacy_spline"
NEW = "tail_spline"
def observation_path() -> Path:
    recovered = RUN / "truth_ionogram_recovered.npz"
    return recovered if recovered.is_file() else RUN / "truth_ionogram.npz"


def layer_from_row(row: dict) -> MonotoneSplineF2:
    return MonotoneSplineF2(row["fof2_mhz"], row["hmf2_km"],
                            row["bottomside_scale_km"],
                            np.asarray(row["log_slope_per_km"], dtype=float))


def write_grid(name: str, layer: MonotoneSplineF2) -> None:
    altitude = np.arange(60.0, 900.1, 20.0)
    latitude = np.arange(-80.0, -27.9, 2.0)
    longitude = np.arange(-24.0, 40.1, 4.0)
    profile = layer.density_cm3(altitude)
    grid = np.broadcast_to(profile, (len(latitude), len(longitude), len(altitude))).copy()
    np.savez_compressed(RUN / f"candidate_{name}_density.npz",
                        latitudes_deg=latitude, longitudes_deg=longitude,
                        altitudes_km=altitude, electron_density_cm3=grid,
                        model=np.array(f"{name} retrieved from ionogram and stated priors"),
                        time_utc=np.array("time withheld during retrieval"))


def fit() -> None:
    path = observation_path()
    observed = Ionogram.read(path)
    measurement_path = RUN / "local_density_observation.json"
    if measurement_path.is_file():
        measurement = json.loads(measurement_path.read_text())
        density = float(measurement["electron_density_cm3"])
        altitude = float(measurement["spacecraft_altitude_km"])
        if not np.isfinite(density) or density <= 0 or altitude != 800.0:
            raise ValueError("Expected a positive local density at the 800 km spacecraft")
        plasma_frequency = PLASMA_MHZ_PER_SQRT_CM3 * np.sqrt(density)
        x_nose, o_nose = observed.nose(-1), observed.nose(1)
        if x_nose is None or o_nose is None:
            raise ValueError("Both O and X returns are needed")
        offsets = ((0.45, 0.35, 0.25) if x_nose - o_nose >= 0.55
                   else (0.25, 0.15, 0.05))
        proposals = (("local_o_tail", True, None),
                     ("local_o_free", False, None),
                     *((f"local_x_{label}", True,
                        max(o_nose + 0.02, x_nose - offset))
                       for label, offset in zip(("low", "mid", "high"), offsets)))
        candidates = {}
        for name, regularize_tail, anchor in proposals:
            result = fit_spline_ionogram(
                observed.records, observed.frequencies, altitude,
                regularize_tail=regularize_tail, fof2_anchor_mhz=anchor,
                spacecraft_plasma_mhz_anchor=plasma_frequency)
            layer = result.layer
            candidates[name] = {
                "fof2_mhz": layer.fof2_mhz,
                "hmf2_km": layer.hmf2_km,
                "bottomside_scale_km": layer.bottomside_scale_km,
                "log_slope_per_km": layer.log_slope_per_km.tolist(),
                "surrogate_range_mae_km": result.range_mae_km[1],
                "tail_regularized": regularize_tail,
                "fof2_anchor_mhz": anchor,
                "spacecraft_density_observation_cm3": density,
                "spacecraft_plasma_frequency_mhz": float(plasma_frequency),
            }
        for name, row in candidates.items():
            write_grid(name, layer_from_row(row))
        document = {
            "observation": path.name,
            "local_density_observation": measurement_path.name,
            "truth_profile_or_peak_parameters_read_during_fit": False,
            "selection_rule": "minimum full-ray O/X accepted-return score among candidates satisfying the local density",
            "candidates": candidates,
        }
        (RUN / "fits.json").write_text(json.dumps(document, indent=2) + "\n")
        print(json.dumps(document, indent=2))
        return
    candidates = {}
    for name, regularize_tail in ((OLD, False), (NEW, True)):
        result = fit_spline_ionogram(observed.records, observed.frequencies,
                                     800.0, regularize_tail=regularize_tail)
        layer = result.layer
        candidates[name] = {
            "fof2_mhz": layer.fof2_mhz,
            "hmf2_km": layer.hmf2_km,
            "bottomside_scale_km": layer.bottomside_scale_km,
            "log_slope_per_km": layer.log_slope_per_km.tolist(),
            "surrogate_range_mae_km": result.range_mae_km[1],
            "tail_regularized": regularize_tail,
        }
    x_nose = observed.nose(-1)
    if x_nose is None:
        raise ValueError("X-mode returns are needed for a mode-resolved cutoff")
    o_nose = observed.nose(1)
    if o_nose is None:
        raise ValueError("O-mode returns are needed for the critical-frequency bound")
    # A genuine O return is a lower bound on foF2. Shift the X-offset grid
    # closer to the X cutoff when the two measured noses are close.
    offsets = ((0.45, 0.35, 0.25) if x_nose - o_nose >= 0.55
               else (0.25, 0.15, 0.05))
    for name, x_offset in zip(("x_anchor_low", "x_anchor_mid", "x_anchor_high"),
                              offsets):
        anchor = max(o_nose + 0.02, x_nose - x_offset)
        result = fit_spline_ionogram(
            observed.records, observed.frequencies, 800.0,
            regularize_tail=True, fof2_anchor_mhz=anchor)
        layer = result.layer
        candidates[name] = {
            "fof2_mhz": layer.fof2_mhz,
            "hmf2_km": layer.hmf2_km,
            "bottomside_scale_km": layer.bottomside_scale_km,
            "log_slope_per_km": layer.log_slope_per_km.tolist(),
            "surrogate_range_mae_km": result.range_mae_km[1],
            "tail_regularized": True,
            "x_nose_mhz": x_nose,
            "assumed_x_cutoff_minus_fof2_mhz": x_offset,
        }
    central_anchor = max(o_nose + 0.02, x_nose - offsets[1])
    for name, upper_frequency in (("tail_freq_04", 0.4),
                                  ("tail_freq_08", 0.8),
                                  ("tail_freq_12", 1.2)):
        result = fit_spline_ionogram(
            observed.records, observed.frequencies, 800.0,
            regularize_tail=True, fof2_anchor_mhz=central_anchor,
            spacecraft_plasma_mhz_anchor=upper_frequency)
        layer = result.layer
        candidates[name] = {
            "fof2_mhz": layer.fof2_mhz,
            "hmf2_km": layer.hmf2_km,
            "bottomside_scale_km": layer.bottomside_scale_km,
            "log_slope_per_km": layer.log_slope_per_km.tolist(),
            "surrogate_range_mae_km": result.range_mae_km[1],
            "tail_regularized": True,
            "x_nose_mhz": x_nose,
            "assumed_x_cutoff_minus_fof2_mhz": offsets[1],
            "spacecraft_plasma_mhz_anchor": upper_frequency,
        }
    for name, row in candidates.items():
        write_grid(name, layer_from_row(row))
    document = {
        "observation": path.name,
        "truth_density_or_peak_parameters_read_during_fit": False,
        "selection_rule": "minimum full-ray O/X accepted-return score",
        "candidates": candidates,
    }
    (RUN / "fits.json").write_text(json.dumps(document, indent=2) + "\n")
    print(json.dumps(document, indent=2))


def height_refine() -> None:
    """Use full-ray range residuals for bounded height steps at fixed local Ne."""
    measurement_path = RUN / "local_density_observation.json"
    if not measurement_path.is_file():
        raise ValueError("Height refinement requires a local density observation")
    document = json.loads((RUN / "fits.json").read_text())
    candidates = document["candidates"]
    observed = Ionogram.read(observation_path())
    first_round = tuple(name for name in candidates if not name.startswith("height_"))
    predictions = {name: Ionogram.read(RUN / f"full_ray_{name}.npz")
                   for name in first_round}
    source_name = min(first_round, key=lambda name: score(observed, predictions[name])["total"])
    source = candidates[source_name]
    residuals = []
    for mode in (1, -1):
        target_ridge = observed.ridge(mode)
        fitted_ridge = predictions[source_name].ridge(mode)
        residuals.extend(fitted_ridge[index] - target_ridge[index]
                         for index in set(target_ridge) & set(fitted_ridge)
                         if observed.frequencies[index] <= 5.0)
    if not residuals:
        raise ValueError("No shared O/X range bins for height refinement")
    median_residual = float(np.median(residuals))
    # A higher reflecting layer shortens the spacecraft-to-reflector path.
    # Thus a positive predicted-minus-observed group range calls for an
    # upward peak-height step in a topside sounder.
    full_step = float(np.clip(median_residual / 2.0, -15.0, 15.0))
    if abs(full_step) < 1.0:
        raise ValueError("Full-ray range residual does not support a height step")
    target_density = float(json.loads(measurement_path.read_text())["electron_density_cm3"])
    for label, factor in (("half", 0.5), ("full", 1.0), ("onehalf", 1.5)):
        name = f"height_{label}"
        new_height = source["hmf2_km"] + factor * full_step
        slopes = np.asarray(source["log_slope_per_km"], dtype=float).copy()

        def local_residual(slope: float) -> float:
            adjusted = slopes.copy()
            adjusted[6] = slope
            layer = MonotoneSplineF2(source["fof2_mhz"], new_height,
                                     source["bottomside_scale_km"], adjusted)
            return float(layer.density_cm3(800.0) - target_density)

        left, right = 0.00002, 0.05
        if local_residual(left) * local_residual(right) >= 0:
            raise ValueError("Cannot preserve local density with one upper slope")
        slopes[6] = brentq(local_residual, left, right, xtol=1e-12)
        row = {**source, "hmf2_km": new_height,
               "log_slope_per_km": slopes.tolist(),
               "source_candidate": source_name,
               "height_step_from_full_ray_km": factor * full_step,
               "spacecraft_density_observation_cm3": target_density}
        candidates[name] = row
        write_grid(name, layer_from_row(row))
    document["height_refinement"] = {
        "source_candidate": source_name,
        "median_2_to_5_mhz_full_ray_range_residual_km": median_residual,
        "full_height_step_km": full_step,
        "measurement_used": measurement_path.name,
        "truth_profile_or_peak_parameters_read": False,
    }
    (RUN / "fits.json").write_text(json.dumps(document, indent=2) + "\n")
    print(json.dumps(document["height_refinement"], indent=2))


def evaluate() -> None:
    observed = Ionogram.read(observation_path())
    gaps = {str(mode): [float(observed.frequencies[index])
                        for group in interior_gaps(observed.records,
                                                   observed.frequencies, mode)
                        for index in group]
            for mode in (1, -1)}
    mode_nose_separation = float(observed.nose(-1) - observed.nose(1))
    observation_qc_warning = (any(gaps.values()) or mode_nose_separation > 0.75)
    candidates = json.loads((RUN / "fits.json").read_text())["candidates"]
    names = tuple(candidates)
    def prediction_path(name: str) -> Path:
        recovered = RUN / f"full_ray_{name}_recovered.npz"
        return recovered if recovered.is_file() else RUN / f"full_ray_{name}.npz"
    predictions = {name: Ionogram.read(prediction_path(name)) for name in names}
    scores = {name: score(observed, ionogram)
              for name, ionogram in predictions.items()}
    winner = min(scores, key=lambda name: scores[name]["total"])
    gap_status = {}
    for name in names:
        with np.load(prediction_path(name), allow_pickle=False) as data:
            details = (json.loads(str(data["gap_recovery_json"]))
                       if "gap_recovery_json" in data else {})
        gap_status[name] = details.get("unrecovered_bins", [])
    (RUN / "selection.json").write_text(json.dumps({
        "observation": observation_path().name,
        "selection_used_truth_density_or_peak_parameters": False,
        "selected": winner,
        "scores": scores,
        "candidate_unrecovered_interior_bins": gap_status,
    }, indent=2) + "\n")
    # The truth files are first opened after the score-only selection above.
    with np.load(RUN / "truth_profile_1km.npz", allow_pickle=False) as data:
        altitude = np.asarray(data["altitudes_km"], dtype=float)
        truth = np.asarray(data["electron_density_cm3"], dtype=float)
    metadata = json.loads((RUN / "provenance.json").read_text())
    truth_fo = metadata["model_fof2_mhz"]
    truth_hm = metadata["model_hmf2_km"]
    upper = (altitude >= truth_hm) & (altitude <= 600.0)
    high_topside = (altitude >= 600.0) & (altitude <= 800.0)
    diagnostics = {}
    for name, row in candidates.items():
        retrieved = layer_from_row(row).density_cm3(altitude)
        diagnostics[name] = {
            "fof2_error_mhz": row["fof2_mhz"] - truth_fo,
            "hmf2_error_km": row["hmf2_km"] - truth_hm,
            "peak_density_error_percent": 100 * ((row["fof2_mhz"] / truth_fo) ** 2 - 1),
            "topside_peak_normalized_rms_percent": float(
                100 * np.sqrt(np.mean((retrieved[upper] - truth[upper]) ** 2))
                / truth.max()),
            "high_topside_relative_rms_percent": float(
                100 * np.sqrt(np.mean(((retrieved[high_topside] - truth[high_topside])
                                       / truth[high_topside]) ** 2))),
            "density_at_800km_cm3": float(retrieved[np.argmin(abs(altitude - 800))]),
            "density_at_800km_error_percent": float(
                100 * (retrieved[np.argmin(abs(altitude - 800))]
                       / truth[np.argmin(abs(altitude - 800))] - 1.0)),
        }
    report = {
        "independent_model": metadata["model"],
        "source_commit": metadata["source_commit"],
        "selection_used_truth_density": False,
        "observation_internal_gaps_mhz": gaps,
        "x_minus_o_nose_mhz": mode_nose_separation,
        "preliminary_due_to_observation_qc": observation_qc_warning,
        "selected": winner,
        "truth_fof2_mhz": truth_fo,
        "truth_hmf2_km": truth_hm,
        "scores": scores,
        "candidate_unrecovered_interior_bins": gap_status,
        "post_selection_density_diagnostics": diagnostics,
    }
    if (RUN / "local_density_observation.json").is_file():
        measurement = json.loads((RUN / "local_density_observation.json").read_text())
        report["local_density_observation_cm3"] = float(measurement["electron_density_cm3"])
        report["selection_used_local_density_observation"] = True
    (RUN / "evaluation.json").write_text(json.dumps(report, indent=2) + "\n")
    plot(observed, predictions[winner], altitude, truth,
         layer_from_row(candidates[winner]).density_cm3(altitude), winner,
         diagnostics[winner])
    print(json.dumps(report, indent=2))


def plot(observed: Ionogram, predicted: Ionogram, altitude: np.ndarray,
         truth: np.ndarray, retrieved: np.ndarray, winner: str, metrics: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(10.7, 7.6), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(1.1, 0.9))
    top_range = max(2150.0, float(np.max(observed.records[:, 2])) + 75.0,
                    float(np.max(predicted.records[:, 2])) + 75.0)
    for ax, ionogram, title in zip(
            (fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])),
            (observed, predicted), ("NeQuick-G independent truth", winner)):
        for mode, color, label in ((1, "#1261ac", "O"), (-1, "#c52525", "X")):
            rows = ionogram.records[ionogram.records[:, 1] == mode]
            ax.scatter(ionogram.frequencies[rows[:, 0].astype(int)],
                       rows[:, 2], color=color, s=12, alpha=0.8,
                       label=f"{label} ({len(rows)})")
        ax.set(xlim=(2, 10), ylim=(150, top_range), xlabel="Frequency (MHz)",
               ylabel="Two-way group range (km)", title=title)
        ax.grid(alpha=0.2)
        ax.legend(frameon=False)
    ax = fig.add_subplot(grid[1, :])
    ax.plot(truth / 1e5, altitude, color="black", linewidth=2.4,
            label="NeQuick-G truth")
    ax.plot(retrieved / 1e5, altitude, color="#c52525", linewidth=2.0,
            label=winner)
    ax.set(xlabel="Electron density (100,000 cm$^{-3}$)",
           ylabel="Altitude (km)", ylim=(150, 800),
           title=(f"foF2 error {metrics['fof2_error_mhz']:+.3f} MHz; "
                  f"hmF2 error {metrics['hmf2_error_km']:+.1f} km; "
                  f"topside RMS {metrics['topside_peak_normalized_rms_percent']:.1f}%"))
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    fig.savefig(RUN / f"{RUN.name}_validation.png", dpi=190)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("fit", "height-refine", "evaluate"))
    parser.add_argument("--run", type=Path, default=RUN,
                        help="Case directory containing saved O/X truth returns")
    args = parser.parse_args()
    RUN = args.run.resolve()
    {"fit": fit, "height-refine": height_refine, "evaluate": evaluate}[args.action]()
