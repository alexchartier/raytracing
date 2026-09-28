"""Evaluate the ionogram-selected pass retrieval against withheld IRI density."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT / "reports") not in sys.path:
    sys.path.insert(0, str(ROOT / "reports"))
from fit_d_ionogram import Ionogram, score  # noqa: E402
from analyze_d_inverse_fit import _fitted_profile  # noqa: E402

DATA = ROOT / "reports/data"
FIGURES = ROOT / "reports/figures"
TRUTH_DENSITY = DATA / "lat_wave_pass_truth_density.npz"
FIT_DENSITY = DATA / "lat_wave_pass_retrieved_density.npz"
TRUTH_IONOGRAMS = DATA / "lat_wave_pass_ionograms_recovered"
FIT_IONOGRAMS = DATA / "lat_wave_pass_retrieved_ionograms"
FIT = DATA / "lat_wave_pass_initial_fit.json"
PRIOR = DATA / "d_inverse_iri_peak_summary.json"
BACKGROUND = DATA / "d_inverse_iri_pyiri_background.npz"
MANIFEST = DATA / "lat_wave_pass_manifest.json"
OUTPUT = DATA / "lat_wave_pass_retrieval_summary.json"


def density_at_positions(path: Path, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as source:
        altitudes = np.asarray(source["altitudes_km"], dtype=float)
        interpolator = RegularGridInterpolator(
            (source["latitudes_deg"], source["longitudes_deg"]),
            source["electron_density_cm3"], bounds_error=True)
        return altitudes, np.asarray(interpolator(points), dtype=float)


def draw_ionogram(ax, path: Path, title: str) -> None:
    with np.load(path, allow_pickle=False) as source:
        records = np.asarray(source["records"], dtype=float)
        frequencies = np.asarray(source["frequencies_mhz"], dtype=float)
    for mode, label, color in ((1, "O", "#135b9a"), (-1, "X", "#b63836")):
        selected = records[records[:, 1] == mode]
        if len(selected):
            # Every accepted ray is drawn at its measured group range. The
            # frequency sampling is 0.1 MHz and range coordinates resolve 1 km.
            ax.scatter(frequencies[selected[:, 0].astype(int)],
                       np.round(selected[:, 2]), s=7, marker="s",
                       linewidths=0, color=color, alpha=0.78, label=label)
    ax.set(xlim=(2, 10), ylim=(150, 2100), title=title)
    ax.grid(alpha=.12)


def evaluate(truth_ionograms: Path = TRUTH_IONOGRAMS,
             fit_ionograms: Path = FIT_IONOGRAMS,
             fit_density: Path = FIT_DENSITY,
             fit_path: Path = FIT,
             output: Path = OUTPUT,
             figure_prefix: str = "lat_wave_retrieval") -> None:
    fit = json.loads(fit_path.read_text())
    manifest = json.loads(MANIFEST.read_text())
    profiles = manifest["profiles"]
    points = np.array([[p["latitude_deg"], p["longitude_deg"]] for p in profiles])
    latitudes = points[:, 0]
    altitudes, truth = density_at_positions(TRUTH_DENSITY, points)
    fit_altitudes, retrieved = density_at_positions(fit_density, points)
    if not np.array_equal(altitudes, fit_altitudes):
        raise ValueError("Truth and retrieved altitude grids differ")
    mask = (altitudes >= 150) & (altitudes <= 600)
    truth_peak = np.max(truth, axis=1)
    fit_peak = np.max(retrieved, axis=1)
    peak_error_pct = 100 * (fit_peak / truth_peak - 1)
    truth_fof2 = .00898 * np.sqrt(truth_peak)
    fit_fof2 = .00898 * np.sqrt(fit_peak)
    profile_nrmse = np.sqrt(np.mean((retrieved[:, mask] - truth[:, mask]) ** 2, axis=1)) / truth_peak
    bottomside = (altitudes >= 150) & (altitudes <= 240)
    topside = (altitudes >= 260) & (altitudes <= 600)
    prior = json.loads(PRIOR.read_text())["retrieved"]
    unperturbed = np.array([_fitted_profile(BACKGROUND, *point, prior)[1]
                            for point in points])
    unperturbed_peak = np.max(unperturbed, axis=1)
    unperturbed_nrmse = np.sqrt(np.mean((unperturbed[:, mask] - truth[:, mask]) ** 2,
                                          axis=1)) / truth_peak
    ionogram_scores = []
    noses = {"truth": [], "retrieved": []}
    counts = {"truth": [], "retrieved": []}
    common_ridge_mae = []
    for p in profiles:
        index = p["index"]
        observed = Ionogram.read(truth_ionograms / f"ionogram_{index:02d}.npz")
        modeled = Ionogram.read(fit_ionograms / f"ionogram_{index:02d}.npz")
        ionogram_scores.append(score(observed, modeled))
        for name, ionogram in (("truth", observed), ("retrieved", modeled)):
            noses[name].append([ionogram.nose(1), ionogram.nose(-1)])
            counts[name].append([int(np.count_nonzero(ionogram.records[:, 1] == mode))
                                 for mode in (1, -1)])
        pair = []
        for mode in (1, -1):
            a, b = observed.ridge(mode), modeled.ridge(mode)
            common = sorted(set(a) & set(b))
            pair.append(float(np.mean([abs(a[i] - b[i]) for i in common])) if common else np.nan)
        common_ridge_mae.append(pair)
    scores = np.array([s["total"] for s in ionogram_scores])
    ridge_mae = np.array(common_ridge_mae)
    truth_noses = np.array(noses["truth"])
    fit_noses = np.array(noses["retrieved"])
    result = {
        "selection_uses_truth_density": False,
        "truth_density_source": "separately generated Fortran IRI-2016 with imposed wave",
        "retrieved_density_source": "PyIRI background transformed by ionogram-only fit",
        "prior_provenance": "PyIRI profile calibrated against one earlier IRI ionogram at pass center; pass density and wave never used for candidate selection",
        "fit_parameters": {key: fit[key] for key in (
            "global_density_scale", "linear_background_fraction_end_to_end_halfspan",
            "wave_amplitude_fraction_at_local_f2_peak", "wave_wavelength_km",
            "wave_phase_rad_at_pass_center", "wave_vertical_sigma_km_assumed")},
        "truth_wave_parameters_post_selection": manifest["wave"],
        "peak_density_median_absolute_error_percent": float(np.median(abs(peak_error_pct))),
        "peak_density_mean_absolute_error_percent": float(np.mean(abs(peak_error_pct))),
        "peak_density_signed_error_percent_by_profile": peak_error_pct.tolist(),
        "fof2_mean_absolute_error_mhz": float(np.mean(abs(fit_fof2 - truth_fof2))),
        "fof2_error_mhz_by_profile": (fit_fof2 - truth_fof2).tolist(),
        "density_profile_mean_nrmse_150_to_600_km": float(np.mean(profile_nrmse)),
        "density_profile_nrmse_by_profile": profile_nrmse.tolist(),
        "bottomside_mean_nrmse_150_to_240_km": float(np.mean(np.sqrt(np.mean(
            (retrieved[:, bottomside] - truth[:, bottomside]) ** 2,
            axis=1)) / truth_peak)),
        "topside_mean_nrmse_260_to_600_km": float(np.mean(np.sqrt(np.mean(
            (retrieved[:, topside] - truth[:, topside]) ** 2,
            axis=1)) / truth_peak)),
        "unperturbed_prior_peak_density_mean_absolute_error_percent": float(
            np.mean(abs(100 * (unperturbed_peak / truth_peak - 1)))),
        "unperturbed_prior_fof2_mean_absolute_error_mhz": float(np.mean(abs(
            .00898 * np.sqrt(unperturbed_peak) - truth_fof2))),
        "unperturbed_prior_density_profile_mean_nrmse_150_to_600_km": float(
            np.mean(unperturbed_nrmse)),
        "ionogram_mean_score": float(np.mean(scores)),
        "ionogram_score_by_profile": scores.tolist(),
        "common_frequency_median_range_mae_km_by_mode":
            np.nanmean(ridge_mae, axis=0).tolist(),
        "nose_mae_mhz_by_mode": np.mean(abs(fit_noses - truth_noses), axis=0).tolist(),
        "noses_mhz_by_profile": noses,
        "accepted_return_counts_by_profile": counts,
        "truth_accepted_returns": int(np.sum(counts["truth"])),
        "retrieved_accepted_returns": int(np.sum(counts["retrieved"])),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")
    FIGURES.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.3), constrained_layout=True,
                             sharex=True, sharey=True)
    mesh_lat, mesh_alt = np.meshgrid(latitudes, altitudes)
    levels = np.linspace(0, max(np.max(truth[:, mask]), np.max(retrieved[:, mask])) / 1e5, 23)
    for ax, values, title in ((axes[0], truth, "Independent truth: IRI-2016 + wave"),
                              (axes[1], retrieved, "Retrieved: PyIRI + fitted wave")):
        image = ax.contourf(mesh_lat, mesh_alt, values.T / 1e5,
                            levels=levels, cmap="viridis", extend="max")
        ax.set(title=title, xlabel="Latitude (°)", ylim=(150, 600))
    fig.colorbar(image, ax=axes[:2], label="Electron density (100,000 cm$^{-3}$)",
                 shrink=.8, pad=.01)
    delta = 100 * (retrieved - truth) / truth_peak[:, None]
    limit = max(10, np.ceil(np.nanmax(abs(delta[:, mask])) / 5) * 5)
    residual = axes[2].contourf(mesh_lat, mesh_alt, delta.T,
                                 levels=np.linspace(-limit, limit, 21),
                                 cmap="RdBu_r", extend="both")
    fig.colorbar(residual, ax=axes[2], label="Difference (% of local truth peak)",
                 shrink=.8, pad=.01)
    axes[2].set(title="Retrieved minus truth", xlabel="Latitude (°)")
    axes[0].set_ylabel("Altitude (km)")
    fig.savefig(FIGURES / f"{figure_prefix}_density.png", dpi=220)
    plt.close(fig)

    indices = (3, 10, 14, 18)
    fig, axes = plt.subplots(len(indices), 2, figsize=(13, 13), constrained_layout=True,
                             sharex=True, sharey=True)
    for row, index in enumerate(indices):
        draw_ionogram(axes[row, 0], truth_ionograms / f"ionogram_{index:02d}.npz",
                      f"Truth, profile {index:02d} ({latitudes[index-1]:.1f}°)")
        draw_ionogram(axes[row, 1], fit_ionograms / f"ionogram_{index:02d}.npz",
                      f"Retrieved, profile {index:02d}")
        axes[row, 0].set_ylabel("Group range (km)")
    for ax in axes[-1]:
        ax.set_xlabel("Frequency (MHz)")
    axes[0, 1].legend(loc="upper right", frameon=True)
    fig.savefig(FIGURES / f"{figure_prefix}_ionograms.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), constrained_layout=True,
                             sharex=True)
    axes[0].plot(latitudes, truth_fof2, "k-", label="Truth density foF2")
    axes[0].plot(latitudes, fit_fof2, color="#24805a", label="Retrieved density foF2")
    for mode, color, marker, label in ((0, "#135b9a", "o", "O return nose"),
                                       (1, "#b63836", "s", "X return nose")):
        axes[0].scatter(latitudes, truth_noses[:, mode], color=color,
                        marker=marker, s=20, alpha=.65, label=label)
        axes[0].plot(latitudes, fit_noses[:, mode], color=color,
                     linestyle="--", linewidth=1.2)
    axes[0].set(ylabel="Frequency (MHz)", title="Density peak and mode-resolved return cutoffs")
    axes[0].legend(ncol=3, fontsize=8)
    axes[1].plot(latitudes, peak_error_pct, "o-", color="#5e3c99")
    axes[1].axhline(0, color="k", linewidth=.8)
    axes[1].set(xlabel="Latitude (°)", ylabel="Peak density error (%)",
                title="Independent peak-density check after ionogram selection")
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(FIGURES / f"{figure_prefix}_peaks.png", dpi=220)
    plt.close(fig)
    print(json.dumps({key: value for key, value in result.items()
                      if not key.endswith("by_profile") and key not in
                      ("noses_mhz_by_profile", "accepted_return_counts_by_profile")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth-ionograms", type=Path, default=TRUTH_IONOGRAMS)
    parser.add_argument("--fit-ionograms", type=Path, default=FIT_IONOGRAMS)
    parser.add_argument("--fit-density", type=Path, default=FIT_DENSITY)
    parser.add_argument("--fit", type=Path, default=FIT)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--figure-prefix", default="lat_wave_retrieval")
    args = parser.parse_args()
    evaluate(args.truth_ionograms, args.fit_ionograms, args.fit_density,
             args.fit, args.output, args.figure_prefix)
