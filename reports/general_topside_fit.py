"""Fit a global-IRI basis to accepted vertical returns and build ray candidates.

`fit` reads only the ionogram and an independent PyIRI background. The IRI-2016
validation density is opened solely by `evaluate`, after candidate selection.
The 1-D fits are proposals; `evaluate` scores saved full 3-D O/X ray returns.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "reports") not in sys.path:
    sys.path.insert(0, str(ROOT / "reports"))

from python_raytrace.general_topside_inverse import (  # noqa: E402
    ChapmanF2, IRITopsideBasis, build_global_iri2020_basis,
    fit_iri_basis_ionogram,
    fit_vertical_ionogram, PLASMA_MHZ_PER_SQRT_CM3,
)
from ionogram_metrics import Ionogram, score  # noqa: E402

DATA = ROOT / "reports/data"
RUN = DATA / "general_topside_iri_case"
BASIS = DATA / "iri2020_global_topside_basis.npz"
OBSERVED = DATA / "d_inverse_iri_truth_blind.npz"
PRIOR = DATA / "d_inverse_iri_pyiri_background.npz"
TRUTH = DATA / "d_inverse_iri_truth_density.npz"
LAT = -53.67859907285381
LON = 7.7400542554165135
SPACECRAFT_ALT_KM = 800.0


def fit() -> None:
    basis = IRITopsideBasis.read(BASIS)
    with np.load(OBSERVED, allow_pickle=False) as data:
        records = np.asarray(data["records"], dtype=float)
        frequencies = np.asarray(data["frequencies_mhz"], dtype=float)
    RUN.mkdir(parents=True, exist_ok=True)
    result = {"basis_source": basis.source,
              "basis_training_profiles": basis.profile_count,
              "observation": str(OBSERVED.relative_to(ROOT)),
              "truth_density_read_during_fit": False,
              "candidate_background": str(PRIOR.relative_to(ROOT)),
              "candidates": {}}
    for name, use_x in (("o_ridge", False), ("ox_ridge", True),
                        ("chapman_o_ridge", False)):
        if name.startswith("chapman"):
            fitted = fit_vertical_ionogram(records, frequencies,
                                           SPACECRAFT_ALT_KM, fit_x_mode=use_x)
        else:
            fitted = fit_iri_basis_ionogram(records, frequencies,
                                            SPACECRAFT_ALT_KM, basis,
                                            fit_x_mode=use_x)
        layer = fitted.layer
        row = {
            "model": "chapman" if isinstance(layer, ChapmanF2) else "iri_basis",
            "fof2_mhz": layer.fof2_mhz,
            "hmf2_km": layer.hmf2_km,
            "x_nuisance_frequency_offset_mhz": fitted.x_frequency_offset_mhz,
            "x_nuisance_range_offset_km": fitted.x_range_offset_km,
            "surrogate_range_mae_km": fitted.range_mae_km,
            "fitted_frequency_bins": fitted.fitted_frequencies,
            "surrogate_objective": fitted.objective,
        }
        if isinstance(layer, ChapmanF2):
            row.update(bottomside_scale_km=layer.bottomside_scale_km,
                       topside_scale_km=layer.topside_scale_km,
                       tail_curvature=layer.tail_curvature)
        else:
            row["basis_scores"] = layer.scores.tolist()
        result["candidates"][name] = row
        build_grid(name, layer)
    destination = RUN / "fits.json"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


def build_grid(name: str, layer, horizontal_gain: float = 1.0) -> None:
    """Carry PyIRI's horizontal peak variation into a flexible F2 profile."""
    if not 0.0 <= horizontal_gain <= 1.5:
        raise ValueError("horizontal_gain must lie in [0, 1.5]")
    with np.load(PRIOR, allow_pickle=False) as source:
        lat = np.asarray(source["latitudes_deg"], dtype=float)
        lon = np.asarray(source["longitudes_deg"], dtype=float)
        alt = np.asarray(source["altitudes_km"], dtype=float)
        prior = np.asarray(source["pyiri_density_cm3"], dtype=float)
    if np.any(prior <= 0) or not np.all(np.isfinite(prior)):
        raise ValueError("PyIRI background has nonpositive density")
    at_sounder = RegularGridInterpolator((lat, lon), prior)([[LAT, LON]])[0]
    reference_fo, reference_hm = peak_from_samples(alt, at_sounder)
    indices = np.argmax(prior, axis=2)
    if np.any((indices == 0) | (indices == len(alt)-1)):
        raise ValueError("PyIRI peak on altitude-grid boundary")
    before = np.take_along_axis(prior, (indices-1)[..., None], axis=2)[..., 0]
    center = np.take_along_axis(prior, indices[..., None], axis=2)[..., 0]
    after = np.take_along_axis(prior, (indices+1)[..., None], axis=2)[..., 0]
    curvature = before - 2*center + after
    if np.any(curvature >= 0):
        raise ValueError("PyIRI contains a non-concave peak")
    spacing = float(alt[1] - alt[0])
    local_hm = alt[indices] + 0.5 * (before - after) / curvature * spacing
    local_peak = center - (after-before)**2 / (8*curvature)
    grid = np.empty_like(prior)
    for i in range(len(lat)):
        for j in range(len(lon)):
            prior_frequency_ratio = (
                PLASMA_MHZ_PER_SQRT_CM3 * np.sqrt(local_peak[i, j])) / reference_fo
            local_fo = layer.fof2_mhz * prior_frequency_ratio ** horizontal_gain
            local_height = layer.hmf2_km + horizontal_gain * (local_hm[i, j] - reference_hm)
            if isinstance(layer, ChapmanF2):
                profile = ChapmanF2(local_fo, local_height,
                                    layer.bottomside_scale_km,
                                    layer.topside_scale_km,
                                    layer.tail_curvature)
            else:
                profile = layer.basis.profile(local_fo, local_height, layer.scores)
            grid[i, j] = profile.density_cm3(alt)
    target = RUN / f"candidate_{name}_density.npz"
    np.savez_compressed(target, latitudes_deg=lat, longitudes_deg=lon,
                        altitudes_km=alt, electron_density_cm3=grid,
                        model=np.array(f"{name} general topside candidate"),
                        time_utc=np.array("2010-01-01T12:00:00"))


def refine() -> None:
    """Propose a small full-ray search from O/X nose discrepancies alone."""
    fitted = json.loads((RUN / "fits.json").read_text())
    base = fitted["candidates"]["o_ridge"]
    basis = IRITopsideBasis.read(BASIS)
    variants = {"basis_fo0995": (0.995, 0.0, 1.0),
                "basis_fo0990": (0.990, 0.0, 1.0),
                "basis_fo0985": (0.985, 0.0, 1.0),
                "basis_fo0990_hmplus5": (0.990, 5.0, 1.0),
                "basis_gradient_015": (1.0, 0.0, 0.15),
                "basis_gradient_030": (1.0, 0.0, 0.30),
                "basis_gradient_050": (1.0, 0.0, 0.5),
                "basis_gradient_000": (1.0, 0.0, 0.0)}
    proposal = {
        "source_candidate": "o_ridge",
        "reason": "Full-ray O/X noses of 5.4/5.7 MHz versus observed 5.4/5.4 MHz; test peak frequency, height, and horizontal prior-gradient gain separately.",
        "selection_uses_truth_density": False,
        "candidates": {},
    }
    for name, (frequency_factor, height_delta, horizontal_gain) in variants.items():
        profile = basis.profile(base["fof2_mhz"] * frequency_factor,
                                base["hmf2_km"] + height_delta,
                                np.asarray(base["basis_scores"]))
        build_grid(name, profile, horizontal_gain=horizontal_gain)
        proposal["candidates"][name] = {
            "model": "iri_basis", "fof2_mhz": profile.fof2_mhz,
            "hmf2_km": profile.hmf2_km,
            "basis_scores": profile.scores.tolist(),
            "frequency_factor_from_o_ridge": frequency_factor,
            "height_delta_from_o_ridge_km": height_delta,
            "horizontal_gradient_gain": horizontal_gain,
        }
    (RUN / "refinements.json").write_text(json.dumps(proposal, indent=2) + "\n")
    print(json.dumps(proposal, indent=2))


def height_refine() -> None:
    """Use shared O/X range residuals to propose a bounded full-ray height step."""
    fitted = json.loads((RUN / "fits.json").read_text())
    base = fitted["candidates"]["o_ridge"]
    observed = Ionogram.read(OBSERVED)
    basis = IRITopsideBasis.read(BASIS)
    proposal = {
        "observation": str(OBSERVED.relative_to(ROOT)),
        "selection_uses_truth_density": False,
        "geometric_two_way_range_derivative_km_per_km": -2.0,
        "sources": {},
        "candidates": {},
    }
    for source_name, prefix, horizontal_gain in (
            ("basis_gradient_000", "basis_flat", 0.0),
            ("basis_gradient_015", "basis_g015", 0.15)):
        predicted = Ionogram.read(RUN / f"full_ray_{source_name}.npz")
        residuals = []
        for mode in (1, -1):
            obs_ridge, pred_ridge = observed.ridge(mode), predicted.ridge(mode)
            for index in sorted(set(obs_ridge) & set(pred_ridge)):
                frequency = observed.frequencies[index]
                if 2.4 <= frequency <= 5.0:
                    residuals.append(pred_ridge[index] - obs_ridge[index])
        if len(residuals) < 20:
            raise ValueError(f"Too few common O/X frequencies for {source_name}")
        median_range_residual = float(np.median(residuals))
        full_step = float(np.clip(median_range_residual / 2.0, -20.0, 20.0))
        proposal["sources"][source_name] = {
            "median_predicted_minus_observed_group_range_km": median_range_residual,
            "full_height_step_km": full_step,
            "common_frequency_count": len(residuals),
        }
        for suffix, factor in (("hm075", 0.75), ("hm100", 1.0)):
            name = f"{prefix}_{suffix}"
            profile = basis.profile(base["fof2_mhz"],
                                    base["hmf2_km"] + factor * full_step,
                                    np.asarray(base["basis_scores"]))
            build_grid(name, profile, horizontal_gain=horizontal_gain)
            proposal["candidates"][name] = {
                "model": "iri_basis", "fof2_mhz": profile.fof2_mhz,
                "hmf2_km": profile.hmf2_km,
                "basis_scores": profile.scores.tolist(),
                "height_step_factor": factor,
                "horizontal_gradient_gain": horizontal_gain,
                "source_candidate": source_name,
            }
    (RUN / "height_refinements.json").write_text(json.dumps(proposal, indent=2) + "\n")
    print(json.dumps(proposal, indent=2))


def evaluate() -> None:
    """Open withheld density only after full-ray selection has been frozen."""
    fitted = json.loads((RUN / "fits.json").read_text())
    observed = Ionogram.read(OBSERVED)
    candidates = dict(fitted["candidates"])
    for filename in ("refinements.json", "height_refinements.json"):
        path = RUN / filename
        if path.exists():
            candidates.update(json.loads(path.read_text())["candidates"])
    options = {}
    for name in candidates:
        path = RUN / f"full_ray_{name}.npz"
        if path.exists():
            options[name] = score(observed, Ionogram.read(path))
    if not options:
        raise FileNotFoundError("No full-ray candidates have been traced")
    selected = min(options, key=lambda name: options[name]["total"])
    basis_options = [name for name in options if candidates[name]["model"] == "iri_basis"]
    best_basis = min(basis_options, key=lambda name: options[name]["total"])
    with np.load(TRUTH, allow_pickle=False) as source:
        altitude = np.asarray(source["altitudes_km"], dtype=float)
        truth = RegularGridInterpolator(
            (source["latitudes_deg"], source["longitudes_deg"]),
            source["electron_density_cm3"])([[LAT, LON]])[0]
    truth_fo, truth_hm = peak_from_samples(altitude, truth)
    top = (altitude >= truth_hm) & (altitude <= 600.0)
    full = (altitude >= 150.0) & (altitude <= 600.0)
    profiles = {}
    diagnostics = {}
    for name in options:
        with np.load(RUN / f"candidate_{name}_density.npz", allow_pickle=False) as source:
            profile = RegularGridInterpolator(
                (source["latitudes_deg"], source["longitudes_deg"]),
                source["electron_density_cm3"])([[LAT, LON]])[0]
        profiles[name] = profile
        fof2, hmf2 = peak_from_samples(altitude, profile)
        diagnostics[name] = {
            "model": candidates[name]["model"],
            "fof2_error_mhz": fof2 - truth_fo,
            "hmf2_error_km": hmf2 - truth_hm,
            "topside_normalized_rms_percent": float(
                100 * np.sqrt(np.mean((truth[top] - profile[top])**2)) / truth.max()),
            "full_150_to_600_normalized_rms_percent": float(
                100 * np.sqrt(np.mean((truth[full] - profile[full])**2)) / truth.max()),
        }
    retrieved = profiles[selected]
    retrieved_fo, retrieved_hm = peak_from_samples(altitude, retrieved)
    summary = {
        "selection": "minimum saved full-ray O/X accepted-return score",
        "selected": selected, "full_ray_scores": options,
        "best_iri_basis_by_full_ray_score": best_basis,
        "post_selection_density_diagnostics": diagnostics,
        "previous_full_ray_score": score(
            observed, Ionogram.read(DATA / "d_inverse_iri_topside_retrieved.npz")),
        "truth_density_source": "IRI-2016 1.11.1",
        "candidate_basis_source": fitted["basis_source"],
        "truth_fof2_mhz": truth_fo, "retrieved_fof2_mhz": retrieved_fo,
        "fof2_error_mhz": retrieved_fo - truth_fo,
        "truth_hmf2_km": truth_hm, "retrieved_hmf2_km": retrieved_hm,
        "hmf2_error_km": retrieved_hm - truth_hm,
        "sampled_truth_hmf2_km": float(altitude[np.argmax(truth)]),
        "sampled_retrieved_hmf2_km": float(altitude[np.argmax(retrieved)]),
        "topside_normalized_rms_percent": diagnostics[selected]["topside_normalized_rms_percent"],
        "full_150_to_600_normalized_rms_percent": diagnostics[selected]["full_150_to_600_normalized_rms_percent"],
        "note": "Peak metrics use three-point quadratic vertices of the 20-km forward grids, not analytic fit parameters.",
    }
    (RUN / "evaluation.json").write_text(json.dumps(summary, indent=2) + "\n")
    plot_evaluation(observed, Ionogram.read(RUN / f"full_ray_{selected}.npz"),
                    altitude, truth, retrieved, summary,
                    candidates[selected]["model"],
                    RUN / "general_topside_full_ray_validation.png")
    if best_basis != selected:
        basis_summary = dict(summary,
                             fof2_error_mhz=diagnostics[best_basis]["fof2_error_mhz"],
                             hmf2_error_km=diagnostics[best_basis]["hmf2_error_km"],
                             topside_normalized_rms_percent=
                             diagnostics[best_basis]["topside_normalized_rms_percent"])
        plot_evaluation(observed, Ionogram.read(RUN / f"full_ray_{best_basis}.npz"),
                        altitude, truth, profiles[best_basis], basis_summary,
                        "iri_basis", RUN / "iri_basis_full_ray_validation.png")
    print(json.dumps(summary, indent=2))


def peak_from_samples(altitude: np.ndarray, density: np.ndarray) -> tuple[float, float]:
    """Three-point peak vertex, consistent with the prior wave evaluation."""
    index = int(np.argmax(density))
    if index == 0 or index == len(altitude) - 1:
        raise ValueError("Peak on altitude-grid boundary")
    before, center, after = density[index-1:index+2]
    curvature = before - 2 * center + after
    if curvature >= 0:
        raise ValueError("Peak is not locally concave")
    spacing = altitude[index+1] - altitude[index]
    shift = 0.5 * (before - after) / curvature * spacing
    peak = center - (after - before)**2 / (8 * curvature)
    return float(PLASMA_MHZ_PER_SQRT_CM3 * np.sqrt(peak)), float(altitude[index] + shift)


def plot_evaluation(observed: Ionogram, predicted: Ionogram,
                    altitude: np.ndarray, truth: np.ndarray,
                    retrieved: np.ndarray, summary: dict, model: str,
                    output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(11.2, 8.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=(1.1, 0.9))
    axes = (fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]))
    for ax, ionogram, title in zip(axes, (observed, predicted),
                                   ("IRI-2016 truth", "Retrieved full-ray candidate")):
        for mode, color, label in ((1, "#1261ac", "O"), (-1, "#c52525", "X")):
            rows = ionogram.records[ionogram.records[:, 1] == mode]
            ax.scatter(ionogram.frequencies[rows[:, 0].astype(int)],
                       np.round(rows[:, 2]),
                       s=11, color=color, alpha=0.75, label=f"{label} ({len(rows)})")
        ax.set(xlim=(2, 10), ylim=(150, 2150), xlabel="Frequency (MHz)",
               ylabel="Two-way group range (km)", title=title)
        ax.grid(alpha=0.2)
        ax.legend(loc="upper left", frameon=False)
    ax = fig.add_subplot(grid[1, :])
    ax.plot(truth / 1e5, altitude, color="black", linewidth=2.3,
            label="IRI-2016 truth")
    ax.plot(retrieved / 1e5, altitude, color="#aa3d1f", linewidth=2.0,
            label=f"Retrieved ({'IRI-2020 basis' if model == 'iri_basis' else 'Chapman'})")
    ax.set(xlabel="Electron density (100,000 cm$^{-3}$)", ylabel="Altitude (km)",
           ylim=(150, 800),
           title=(f"foF2 error {summary['fof2_error_mhz']:+.3f} MHz; "
                  f"hmF2 error {summary['hmf2_error_km']:+.0f} km; "
                  f"topside RMS {summary['topside_normalized_rms_percent']:.2f}%"))
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    fig.savefig(output, dpi=190)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("build-basis", "fit", "refine",
                                           "height-refine", "evaluate"))
    args = parser.parse_args()
    if args.action == "build-basis":
        basis = build_global_iri2020_basis()
        basis.save(BASIS)
        print(json.dumps({"basis": str(BASIS), "source": basis.source,
                          "profile_count": basis.profile_count}))
    elif args.action == "fit":
        fit()
    elif args.action == "refine":
        refine()
    elif args.action == "height-refine":
        height_refine()
    else:
        evaluate()
