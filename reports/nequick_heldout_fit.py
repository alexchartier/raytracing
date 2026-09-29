"""Fit a second NeQuick-G ionogram using only its saved O/X returns.

The retrieval never imports NeQuick. The independent density and its peak
parameters are opened only by ``evaluate``, after full-ray score selection.
This July profile became a development case when its O/X cutoff difference
was inspected; do not label it a prospective held-out validation of that rule.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "reports"))

from ionogram_metrics import Ionogram, score  # noqa: E402
from python_raytrace.general_topside_inverse import (  # noqa: E402
    MonotoneSplineF2, fit_spline_ionogram,
)

RUN = ROOT / "reports/data/nequick_heldout_2026_case"
OLD = "legacy_spline"
NEW = "tail_spline"
NAMES = (OLD, NEW, "x_anchor_low", "x_anchor_mid", "x_anchor_high")


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
                        model=np.array(f"{name} retrieved from O-mode returns"),
                        time_utc=np.array("month 07, 8 UT"))


def fit() -> None:
    observed = Ionogram.read(RUN / "truth_ionogram.npz")
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
    for name, x_offset in (("x_anchor_low", 0.45),
                           ("x_anchor_mid", 0.35),
                           ("x_anchor_high", 0.25)):
        anchor = x_nose - x_offset
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
    for name, row in candidates.items():
        write_grid(name, layer_from_row(row))
    document = {
        "observation": "truth_ionogram.npz",
        "truth_density_or_peak_parameters_read_during_fit": False,
        "selection_rule": "minimum full-ray O/X accepted-return score",
        "candidates": candidates,
    }
    (RUN / "fits.json").write_text(json.dumps(document, indent=2) + "\n")
    print(json.dumps(document, indent=2))


def evaluate() -> None:
    observed = Ionogram.read(RUN / "truth_ionogram.npz")
    candidates = json.loads((RUN / "fits.json").read_text())["candidates"]
    predictions = {name: Ionogram.read(RUN / f"full_ray_{name}.npz")
                   for name in NAMES}
    scores = {name: score(observed, ionogram)
              for name, ionogram in predictions.items()}
    winner = min(scores, key=lambda name: scores[name]["total"])
    # The truth files are first opened after the score-only selection above.
    with np.load(RUN / "truth_profile_1km.npz", allow_pickle=False) as data:
        altitude = np.asarray(data["altitudes_km"], dtype=float)
        truth = np.asarray(data["electron_density_cm3"], dtype=float)
    metadata = json.loads((RUN / "provenance.json").read_text())
    truth_fo = metadata["model_fof2_mhz"]
    truth_hm = metadata["model_hmf2_km"]
    upper = (altitude >= truth_hm) & (altitude <= 600.0)
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
            "density_at_800km_cm3": float(retrieved[np.argmin(abs(altitude - 800))]),
        }
    report = {
        "independent_model": metadata["model"],
        "source_commit": metadata["source_commit"],
        "selection_used_truth_density": False,
        "selected": winner,
        "truth_fof2_mhz": truth_fo,
        "truth_hmf2_km": truth_hm,
        "scores": scores,
        "post_selection_density_diagnostics": diagnostics,
    }
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
    for ax, ionogram, title in zip(
            (fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])),
            (observed, predicted), ("NeQuick-G held-out truth", winner)):
        for mode, color, label in ((1, "#1261ac", "O"), (-1, "#c52525", "X")):
            rows = ionogram.records[ionogram.records[:, 1] == mode]
            ax.scatter(ionogram.frequencies[rows[:, 0].astype(int)],
                       rows[:, 2], color=color, s=12, alpha=0.8,
                       label=f"{label} ({len(rows)})")
        ax.set(xlim=(2, 10), ylim=(150, 2150), xlabel="Frequency (MHz)",
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
    fig.savefig(RUN / "heldout_nequick_validation.png", dpi=190)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("fit", "evaluate"))
    {"fit": fit, "evaluate": evaluate}[parser.parse_args().action]()
