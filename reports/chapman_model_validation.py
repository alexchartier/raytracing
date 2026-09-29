"""Fit analytic-Chapman O/X truth without reading its generating profile.

`fit` reads only the saved O/X ionogram and precomputed IRI basis. `evaluate`
opens truth density and parameters after full-ray scores select a candidate.
The separate `prepare_chapman_truth.py` generates truth without retrieval code.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "reports"))

from ionogram_metrics import Ionogram, score  # noqa: E402
from python_raytrace.general_topside_inverse import (  # noqa: E402
    ChapmanF2, IRITopsideBasis, fit_iri_basis_ionogram, fit_vertical_ionogram,
)

DATA = ROOT / "reports/data"
OUT = DATA / "chapman_model_validation"
BASIS = DATA / "iri2020_global_topside_basis.npz"
TRUTH = OUT / "truth_density.npz"
OBSERVED = OUT / "truth_ionogram.npz"
LAT = np.arange(-80.0, -28.0 + 0.1, 2.0)
LON = np.arange(-24.0, 40.0 + 0.1, 4.0)
ALT = np.arange(60.0, 900.0 + 0.1, 20.0)
NAMES = ("chapman_o", "chapman_ox", "iri_basis",
         "chapman_o_hm20", "chapman_o_hm40", "chapman_o_refined")


def write_grid(path: Path, lat: np.ndarray, lon: np.ndarray,
               alt: np.ndarray, layer, label: str) -> None:
    profile = layer.density_cm3(alt)
    grid = np.broadcast_to(profile, (len(lat), len(lon), len(alt))).copy()
    np.savez_compressed(path, latitudes_deg=lat, longitudes_deg=lon,
                        altitudes_km=alt, electron_density_cm3=grid,
                        model=np.array(label), time_utc=np.array("2010-01-01T12:00:00"))


def fit() -> None:
    observed = Ionogram.read(OBSERVED)
    basis = IRITopsideBasis.read(BASIS)
    layers = {
        "chapman_o": fit_vertical_ionogram(observed.records, observed.frequencies,
                                            800.0, fit_x_mode=False),
        "chapman_ox": fit_vertical_ionogram(observed.records, observed.frequencies,
                                             800.0, fit_x_mode=True),
        "iri_basis": fit_iri_basis_ionogram(observed.records, observed.frequencies,
                                             800.0, basis, fit_x_mode=True),
    }
    candidate_layers = {name: result.layer for name, result in layers.items()}
    for delta in (20.0, 40.0):
        candidate_layers[f"chapman_o_hm{int(delta)}"] = replace(
            layers["chapman_o"].layer,
            hmf2_km=layers["chapman_o"].layer.hmf2_km + delta)
    candidates = {}
    for name, layer in candidate_layers.items():
        write_grid(OUT / f"candidate_{name}_density.npz", LAT, LON, ALT,
                   layer, f"{name} retrieved from O/X ionogram")
        source = "chapman_o" if name.startswith("chapman_o_hm") else name
        result = layers[source]
        candidates[name] = {
            "fof2_mhz": layer.fof2_mhz,
            "hmf2_km": layer.hmf2_km,
            "source_fit": source,
            "height_step_from_source_km": layer.hmf2_km - result.layer.hmf2_km,
            "surrogate_range_mae_km": result.range_mae_km,
            "surrogate_objective": result.objective,
            "x_nuisance_frequency_offset_mhz": result.x_frequency_offset_mhz,
            "x_nuisance_range_offset_km": result.x_range_offset_km,
            **({"topside_scale_km": layer.topside_scale_km,
                "tail_curvature": layer.tail_curvature}
               if isinstance(layer, ChapmanF2) else
               {"basis_scores": layer.scores.tolist()}),
        }
    (OUT / "fits.json").write_text(json.dumps({
        "fitted_from": str(OBSERVED.relative_to(ROOT)),
        "truth_density_read_during_fit": False,
        "candidates": candidates,
    }, indent=2) + "\n")
    print(json.dumps(candidates, indent=2))


def refine() -> None:
    """Convert a robust full-ray group-range residual into a height proposal."""
    observed = Ionogram.read(OBSERVED)
    predicted = Ionogram.read(OUT / "full_ray_chapman_o.npz")
    first = min(observed.frequencies[observed.records[:, 0].astype(int)])
    upper = observed.nose(1) - 0.5
    lower = max(2.4, first + 0.5)
    residuals = []
    for mode in (1, -1):
        actual, forecast = observed.ridge(mode), predicted.ridge(mode)
        for index in sorted(set(actual) & set(forecast)):
            if lower <= observed.frequencies[index] <= upper:
                residuals.append(forecast[index] - actual[index])
    if len(residuals) < 20:
        raise ValueError("Too few shared O/X bins for a height refinement")
    median = float(np.median(residuals))
    height_step = float(np.clip(median / 2.0, -60.0, 60.0))
    fitted = json.loads((OUT / "fits.json").read_text())
    base = fitted["candidates"]["chapman_o"]
    profile = ChapmanF2(base["fof2_mhz"], base["hmf2_km"] + height_step,
                        55.0, base["topside_scale_km"], base["tail_curvature"])
    write_grid(OUT / "candidate_chapman_o_refined_density.npz", LAT, LON, ALT,
               profile, "Chapman candidate refined from full-ray O/X ranges")
    fitted["candidates"]["chapman_o_refined"] = {
        **base, "hmf2_km": profile.hmf2_km,
        "source_fit": "chapman_o",
        "height_step_from_source_km": height_step,
        "median_predicted_minus_observed_range_km": median,
        "shared_frequency_count": len(residuals),
        "refinement_frequency_window_mhz": [lower, upper],
        "selection_uses_truth_density": False,
    }
    (OUT / "fits.json").write_text(json.dumps(fitted, indent=2) + "\n")
    print(json.dumps(fitted["candidates"]["chapman_o_refined"], indent=2))


def evaluate() -> None:
    observed = Ionogram.read(OBSERVED)
    scores = {name: score(observed, Ionogram.read(OUT / f"full_ray_{name}.npz"))
              for name in NAMES}
    selected = min(scores, key=lambda name: scores[name]["total"])
    with np.load(TRUTH, allow_pickle=False) as source:
        alt = np.asarray(source["altitudes_km"], dtype=float)
        truth = np.asarray(source["electron_density_cm3"][0, 0], dtype=float)
    truth_parameters = json.loads((OUT / "truth_parameters.json").read_text())
    fitted = json.loads((OUT / "fits.json").read_text())["candidates"]
    truth_fo = float(truth_parameters["fof2_mhz"])
    truth_hm = float(truth_parameters["hmf2_km"])
    top = (alt >= truth_hm) & (alt <= 600.0)
    diagnostics = {}
    profiles = {}
    for name in NAMES:
        with np.load(OUT / f"candidate_{name}_density.npz", allow_pickle=False) as source:
            profile = np.asarray(source["electron_density_cm3"][0, 0], dtype=float)
        profiles[name] = profile
        fof2 = float(fitted[name]["fof2_mhz"])
        hmf2 = float(fitted[name]["hmf2_km"])
        diagnostics[name] = {
            "fof2_error_mhz": fof2 - truth_fo,
            "peak_density_error_percent": 100 * ((fof2 / truth_fo)**2 - 1),
            "hmf2_error_km": hmf2 - truth_hm,
            "topside_normalized_rms_percent": float(
                100 * np.sqrt(np.mean((profile[top] - truth[top])**2)) / np.max(truth)),
        }
    summary = {
        "truth_source": "analytic generalized Chapman; uniform 3-D field",
        "selection": "minimum full-ray O/X accepted-return score",
        "selected": selected,
        "truth_fof2_mhz": truth_fo,
        "truth_hmf2_km": truth_hm,
        "peak_metric_definition": "continuous analytic truth and fitted profile parameters; 20-km grid is used only for the RMS diagnostic",
        "scores": scores,
        "post_selection_density_diagnostics": diagnostics,
    }
    (OUT / "evaluation.json").write_text(json.dumps(summary, indent=2) + "\n")
    plot(observed, Ionogram.read(OUT / f"full_ray_{selected}.npz"),
         alt, truth, profiles, selected, summary)
    print(json.dumps(summary, indent=2))


def plot(observed: Ionogram, predicted: Ionogram, alt: np.ndarray,
         truth: np.ndarray, profiles: dict, selected: str, summary: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(10.7, 7.6), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(1.1, 0.9))
    for ax, ionogram, title in zip(
            (fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])),
            (observed, predicted), ("Chapman truth", "Selected Chapman retrieval")):
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
    ax.plot(truth / 1e5, alt, color="black", linewidth=2.4, label="Chapman truth")
    comparison = [(selected, "#c52525", "-")]
    if selected != "iri_basis":
        comparison.append(("iri_basis", "#1261ac", "--"))
    for name, color, style in comparison:
        label = "Selected Chapman" if name == selected else "Forced IRI basis"
        ax.plot(profiles[name] / 1e5, alt, color=color, linestyle=style,
                linewidth=2.0, label=label)
    ax.set(xlabel="Electron density (100,000 cm$^{-3}$)",
           ylabel="Altitude (km)", ylim=(150, 800),
           title=("foF2 error "
                  f"{summary['post_selection_density_diagnostics'][selected]['fof2_error_mhz']:+.3f} MHz; "
                  "hmF2 error "
                  f"{summary['post_selection_density_diagnostics'][selected]['hmf2_error_km']:+.1f} km"))
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    fig.savefig(OUT / "chapman_full_ray_validation.png", dpi=190)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("fit", "refine", "evaluate"))
    action = parser.parse_args().action
    {"fit": fit, "refine": refine, "evaluate": evaluate}[action]()
