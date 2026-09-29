"""Fit saved NeQuick-G O/X returns without loading the NeQuick-G model or density.

The input families are an IRI-2020 basis, a Chapman layer, and an unconstrained
monotone log-density spline. Only `evaluate` reads the independent truth.
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
    ChapmanF2, IRITopsideBasis, MonotoneSplineF2, fit_iri_basis_ionogram,
    fit_spline_ionogram, fit_vertical_ionogram,
)

DATA = ROOT / "reports/data"
OUT = DATA / "nequick_independent_case"
OBSERVED = OUT / "truth_ionogram.npz"
BASIS = DATA / "iri2020_global_topside_basis.npz"
LAT = np.arange(-80.0, -28.0 + 0.1, 2.0)
LON = np.arange(-24.0, 40.0 + 0.1, 4.0)
ALT = np.arange(60.0, 900.0 + 0.1, 20.0)
NAMES = ("iri_basis", "chapman", "flexible", "flexible_refined")


def write_grid(name: str, layer) -> None:
    profile = layer.density_cm3(ALT)
    density = np.broadcast_to(profile, (len(LAT), len(LON), len(ALT))).copy()
    np.savez_compressed(OUT / f"candidate_{name}_density.npz",
                        latitudes_deg=LAT, longitudes_deg=LON,
                        altitudes_km=ALT, electron_density_cm3=density,
                        model=np.array(f"{name} candidate from NeQuick-G ionogram"),
                        time_utc=np.array("month 01, 12 UT"))


def fit() -> None:
    observed = Ionogram.read(OBSERVED)
    basis = IRITopsideBasis.read(BASIS)
    fitted = {
        "iri_basis": fit_iri_basis_ionogram(observed.records,
                                             observed.frequencies, 800.0,
                                             basis, fit_x_mode=False),
        "chapman": fit_vertical_ionogram(observed.records,
                                          observed.frequencies, 800.0,
                                          fit_x_mode=False),
        "flexible": fit_spline_ionogram(observed.records,
                                        observed.frequencies, 800.0),
    }
    rows = {}
    for name, result in fitted.items():
        layer = result.layer
        write_grid(name, layer)
        parameters = {
            "fof2_mhz": layer.fof2_mhz,
            "hmf2_km": layer.hmf2_km,
            "bottomside_scale_km": layer.bottomside_scale_km,
            "surrogate_range_mae_km": result.range_mae_km,
            "surrogate_objective": result.objective,
        }
        if isinstance(layer, ChapmanF2):
            parameters.update(topside_scale_km=layer.topside_scale_km,
                              tail_curvature=layer.tail_curvature)
        elif isinstance(layer, MonotoneSplineF2):
            parameters["log_slope_per_km"] = layer.log_slope_per_km.tolist()
        else:
            parameters["basis_scores"] = layer.scores.tolist()
        rows[name] = parameters
    (OUT / "fits.json").write_text(json.dumps({
        "observation": str(OBSERVED.relative_to(ROOT)),
        "truth_density_or_parameters_read_during_fit": False,
        "candidates": rows,
    }, indent=2) + "\n")
    print(json.dumps(rows, indent=2))


def refine() -> None:
    observed = Ionogram.read(OBSERVED)
    predicted = Ionogram.read(OUT / "full_ray_flexible.npz")
    upper = observed.nose(1) - 0.4
    lower = max(2.4, min(observed.frequencies[observed.records[:, 0].astype(int)]) + 0.4)
    residuals = []
    for mode in (1, -1):
        actual, forecast = observed.ridge(mode), predicted.ridge(mode)
        for index in sorted(set(actual) & set(forecast)):
            if lower <= observed.frequencies[index] <= upper:
                residuals.append(forecast[index] - actual[index])
    if len(residuals) < 15:
        raise ValueError("Too few common O/X bins for a height refinement")
    median = float(np.median(residuals))
    step = float(np.clip(median / 2.0, -40.0, 40.0))
    fitted = json.loads((OUT / "fits.json").read_text())
    source = fitted["candidates"]["flexible"]
    layer = MonotoneSplineF2(source["fof2_mhz"], source["hmf2_km"] + step,
                             source["bottomside_scale_km"],
                             np.asarray(source["log_slope_per_km"], dtype=float))
    write_grid("flexible_refined", layer)
    fitted["candidates"]["flexible_refined"] = {
        **source, "hmf2_km": layer.hmf2_km,
        "height_step_from_full_ray_residual_km": step,
        "median_predicted_minus_observed_range_km": median,
        "common_frequency_count": len(residuals),
        "frequency_window_mhz": [lower, upper],
        "truth_density_or_parameters_read_during_refinement": False,
    }
    (OUT / "fits.json").write_text(json.dumps(fitted, indent=2) + "\n")
    print(json.dumps(fitted["candidates"]["flexible_refined"], indent=2))


def reconstruct(name: str, parameters: dict, basis: IRITopsideBasis):
    if name == "iri_basis":
        return basis.profile(parameters["fof2_mhz"], parameters["hmf2_km"],
                             np.asarray(parameters["basis_scores"]),
                             parameters["bottomside_scale_km"])
    if name == "chapman":
        return ChapmanF2(parameters["fof2_mhz"], parameters["hmf2_km"],
                         parameters["bottomside_scale_km"],
                         parameters["topside_scale_km"],
                         parameters["tail_curvature"])
    return MonotoneSplineF2(parameters["fof2_mhz"], parameters["hmf2_km"],
                            parameters["bottomside_scale_km"],
                            np.asarray(parameters["log_slope_per_km"]))


def evaluate() -> None:
    if (OUT / "truth_ionogram_x_gap_recovered.npz").is_file():
        raise RuntimeError(
            "The original adaptive truth has a recovered X gap. Re-run candidate "
            "ionograms with the same gap recovery before scoring them against it."
        )
    observed = Ionogram.read(OBSERVED)
    available = [name for name in NAMES if (OUT / f"full_ray_{name}.npz").is_file()]
    if not {"iri_basis", "flexible", "flexible_refined"}.issubset(available):
        raise FileNotFoundError("The required IRI and flexible full-ray checks are incomplete")
    scores = {name: score(observed, Ionogram.read(OUT / f"full_ray_{name}.npz"))
              for name in available}
    selected = min(scores, key=lambda name: scores[name]["total"])
    # All truth reads follow full-ray selection.
    with np.load(OUT / "truth_profile_1km.npz", allow_pickle=False) as source:
        altitude = np.asarray(source["altitudes_km"], dtype=float)
        truth = np.asarray(source["electron_density_cm3"], dtype=float)
    truth_metadata = json.loads((OUT / "provenance.json").read_text())
    truth_fo = float(truth_metadata["model_fof2_mhz"])
    truth_hm = float(truth_metadata["model_hmf2_km"])
    fitted = json.loads((OUT / "fits.json").read_text())["candidates"]
    basis = IRITopsideBasis.read(BASIS)
    top = (altitude >= truth_hm) & (altitude <= 600.0)
    diagnostics = {}
    profiles = {}
    for name in NAMES:
        profile = reconstruct(name, fitted[name], basis).density_cm3(altitude)
        profiles[name] = profile
        fo = fitted[name]["fof2_mhz"]
        diagnostics[name] = {
            "fof2_error_mhz": fo - truth_fo,
            "hmf2_error_km": fitted[name]["hmf2_km"] - truth_hm,
            "peak_density_error_percent": 100 * ((fo / truth_fo)**2 - 1),
            "topside_normalized_rms_percent": float(
                100 * np.sqrt(np.mean((profile[top] - truth[top])**2)) / truth.max()),
        }
    summary = {
        "truth_model": truth_metadata["model"],
        "source_commit": truth_metadata["source_commit"],
        "selection": "minimum completed full-ray O/X accepted-return score",
        "selected": selected,
        "unscored_candidates": sorted(set(NAMES) - set(available)),
        "truth_fof2_mhz": truth_fo,
        "truth_hmf2_km": truth_hm,
        "scores": scores,
        "post_selection_density_diagnostics": diagnostics,
    }
    (OUT / "evaluation.json").write_text(json.dumps(summary, indent=2) + "\n")
    plot(observed, Ionogram.read(OUT / f"full_ray_{selected}.npz"),
         altitude, truth, profiles[selected], selected, diagnostics[selected])
    print(json.dumps(summary, indent=2))


def plot(observed: Ionogram, predicted: Ionogram, altitude: np.ndarray,
         truth: np.ndarray, retrieved: np.ndarray, name: str, diagnostics: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(10.7, 7.6), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(1.1, 0.9))
    display_name = "Free-form spline retrieval" if name.startswith("flexible") else name
    for ax, ionogram, title in zip(
            (fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])),
            (observed, predicted), ("NeQuick-G truth", display_name)):
        for mode, color, label in ((1, "#1261ac", "O"), (-1, "#c52525", "X")):
            rows = ionogram.records[ionogram.records[:, 1] == mode]
            ax.scatter(ionogram.frequencies[rows[:, 0].astype(int)],
                       np.round(rows[:, 2]), color=color, s=12, alpha=0.8,
                       label=f"{label} ({len(rows)})")
        ax.set(xlim=(2, 10), ylim=(150, 2150), xlabel="Frequency (MHz)",
               ylabel="Two-way group range (km)", title=title)
        ax.grid(alpha=0.2)
        ax.legend(frameon=False)
    ax = fig.add_subplot(grid[1, :])
    ax.plot(truth / 1e5, altitude, color="black", linewidth=2.4,
            label="NeQuick-G truth")
    ax.plot(retrieved / 1e5, altitude, color="#c52525", linewidth=2.0,
            label=display_name)
    ax.set(xlabel="Electron density (100,000 cm$^{-3}$)",
           ylabel="Altitude (km)", ylim=(150, 800),
           title=(f"foF2 error {diagnostics['fof2_error_mhz']:+.3f} MHz; "
                  f"hmF2 error {diagnostics['hmf2_error_km']:+.1f} km; "
                  f"topside RMS {diagnostics['topside_normalized_rms_percent']:.1f}%"))
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    fig.savefig(OUT / "nequick_full_ray_validation.png", dpi=190)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("fit", "refine", "evaluate"))
    action = parser.parse_args().action
    {"fit": fit, "refine": refine, "evaluate": evaluate}[action]()
