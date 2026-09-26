"""Summarize a completed D ionogram fit and, optionally, its density profile.

The profile comparison uses a separately generated truth grid only after the
ionogram search has finished. It never enters candidate scoring or selection.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

from fit_d_ionogram import Ionogram


def _profile(grid_path: Path, field: str, latitude: float, longitude: float) -> tuple[np.ndarray, np.ndarray]:
    with np.load(grid_path, allow_pickle=False) as data:
        altitudes = np.asarray(data["altitudes_km"], dtype=float)
        interpolator = RegularGridInterpolator(
            (data["latitudes_deg"], data["longitudes_deg"]),
            data[field], bounds_error=True,
        )
        profile = np.asarray(interpolator([[latitude, longitude]])[0], dtype=float)
    return altitudes, profile


def _fitted_profile(background_path: Path, latitude: float, longitude: float,
                    candidate: dict) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate the candidate's transformed grid at the sounder location."""
    with np.load(background_path, allow_pickle=False) as background:
        altitudes = np.asarray(background["altitudes_km"], dtype=float)
        lats = np.asarray(background["latitudes_deg"], dtype=float)
        lons = np.asarray(background["longitudes_deg"], dtype=float)
        grid = np.asarray(background["pyiri_density_cm3"], dtype=float)
    i = int(np.searchsorted(lats, latitude) - 1)
    j = int(np.searchsorted(lons, longitude) - 1)
    if not 0 <= i < len(lats) - 1 or not 0 <= j < len(lons) - 1:
        raise ValueError("Sounder lies outside the background grid")
    u = (latitude - lats[i]) / (lats[i + 1] - lats[i])
    v = (longitude - lons[j]) / (lons[j + 1] - lons[j])
    width = float(candidate.get("f2_width_scale", 1.0))
    shift = float(candidate["hmf2_shift_km"])
    result = np.zeros_like(altitudes)
    for di, lat_weight in ((0, 1-u), (1, u)):
        for dj, lon_weight in ((0, 1-v), (1, v)):
            base = grid[i + di, j + dj]
            peak_altitude = float(altitudes[np.argmax(base)])
            mapped = peak_altitude + (altitudes - peak_altitude) / width
            stretched = np.interp(mapped, altitudes, base, left=base[0], right=base[-1])
            shifted = np.interp(altitudes - shift, altitudes, stretched,
                                left=stretched[0], right=stretched[-1])
            result += lat_weight * lon_weight * shifted
    return altitudes, float(candidate["density_scale"]) * result


def summarize(state_path: Path, output_prefix: Path, background_path: Path | None,
              truth_density_path: Path | None, figure_output: Path | None = None) -> None:
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    state = json.loads(state_path.read_text())
    evaluations = state["evaluations"]
    best = min(evaluations, key=lambda row: row["score"]["total"])
    observed = Ionogram.read(Path(state["observed"]))
    retrieved = Ionogram.read(Path(best["path"]))
    summary = {
        "candidate_count": len(evaluations),
        "best_score_by_population": [
            float(min(row["score"]["total"] for row in evaluations[:stop]))
            for stop in (40, 80, 120, 160) if stop <= len(evaluations)
        ],
        "retrieved": {key: best[key] for key in ("density_scale", "hmf2_shift_km", "score")},
        "returns": {
            label: {str(mode): int(np.count_nonzero(ionogram.records[:, 1] == mode))
                    for mode in (1, -1)}
            for label, ionogram in (("truth", observed), ("retrieved", retrieved))
        },
        "nose_mhz": {
            label: {str(mode): ionogram.nose(mode) for mode in (1, -1)}
            for label, ionogram in (("truth", observed), ("retrieved", retrieved))
        },
        "common_frequency_median_range_mae_km": {},
    }
    if "f2_width_scale" in best:
        summary["retrieved"]["f2_width_scale"] = best["f2_width_scale"]
    for mode in (1, -1):
        truth_ridge, fit_ridge = observed.ridge(mode), retrieved.ridge(mode)
        common = sorted(set(truth_ridge) & set(fit_ridge))
        summary["common_frequency_median_range_mae_km"][str(mode)] = (
            float(np.mean([abs(truth_ridge[index] - fit_ridge[index]) for index in common]))
            if common else None
        )
    if (background_path is None) != (truth_density_path is None):
        raise ValueError("Supply both the background and independent truth grids")
    if background_path is not None and truth_density_path is not None:
        with np.load(background_path, allow_pickle=False) as background:
            latitude = float(background["tx_lat_deg"])
            longitude = float(background["tx_lon_deg"])
        width = float(best.get("f2_width_scale", 1.0))
        altitudes, fit = _fitted_profile(background_path, latitude, longitude, best)
        truth_altitudes, truth = _profile(truth_density_path, "electron_density_cm3", latitude, longitude)
        if not np.array_equal(altitudes, truth_altitudes):
            raise ValueError("Density altitude grids differ")
        mask = (altitudes >= 150) & (altitudes <= 600)
        peak = float(np.max(truth[mask]))
        peak_altitude = float(altitudes[np.argmax(truth)])
        bottomside = (altitudes >= 150) & (altitudes < peak_altitude)
        topside = (altitudes >= peak_altitude) & (altitudes <= 600)
        def profile_error(row: dict) -> float:
            candidate_profile = _fitted_profile(background_path, latitude, longitude, row)[1]
            return float(np.sqrt(np.mean((candidate_profile[mask] - truth[mask]) ** 2)) / peak)

        oracle = min(
            ((profile_error(row), row) for row in evaluations),
            key=lambda pair: pair[0],
        )
        oracle_profile = _fitted_profile(background_path, latitude, longitude, oracle[1])[1]
        summary["density_at_sounder"] = {
            "latitude_deg": latitude, "longitude_deg": longitude,
            "truth_peak_cm3": float(np.max(truth)),
            "truth_peak_altitude_km": peak_altitude,
            "retrieved_peak_cm3": float(np.max(fit)),
            "retrieved_peak_altitude_km": float(altitudes[np.argmax(fit)]),
            "rms_error_150_to_600_km_relative_to_truth_peak": float(
                np.sqrt(np.mean((truth[mask] - fit[mask]) ** 2)) / peak),
            "median_absolute_percentage_error_150_to_600_km": float(
                np.median(np.abs(truth[mask] - fit[mask]) / truth[mask]) * 100),
            "rms_error_peak_to_600_km_relative_to_truth_peak": float(
                np.sqrt(np.mean((truth[topside] - fit[topside]) ** 2)) / peak),
            "median_absolute_percentage_error_peak_to_600_km": float(
                np.median(np.abs(truth[topside] - fit[topside]) / truth[topside]) * 100),
            "rms_error_150_km_to_below_peak_relative_to_truth_peak": float(
                np.sqrt(np.mean((truth[bottomside] - fit[bottomside]) ** 2)) / peak),
        }
        summary["profile_oracle_diagnostic"] = {
            "uses_truth_density_for_selection": True,
            "candidate": {key: oracle[1].get(key, 1.0) for key in
                          ("density_scale", "hmf2_shift_km", "f2_width_scale")},
            "ionogram_score": oracle[1]["score"]["total"],
            "rms_error_150_to_600_km_relative_to_truth_peak": oracle[0],
        }
        fig, ax = plt.subplots(figsize=(7, 6), constrained_layout=True)
        ax.plot(truth / 1e5, altitudes, color="#9a4522", linewidth=2.5, label="IRI-2016 truth")
        ax.plot(fit / 1e5, altitudes, color="#253a5e", linewidth=2.5,
                label=("Retrieved PyIRI scale + height + width" if width != 1.0
                       else "Retrieved PyIRI scale + height"))
        if oracle[1]["path"] != best["path"]:
            ax.plot(oracle_profile / 1e5, altitudes, color="#4d7a55", linewidth=2,
                    linestyle="--", label="Best profile among candidates (truth selected)")
        ax.axhline(peak_altitude, color="#777777", linestyle=":", linewidth=1)
        ax.set(xlabel="Electron density (100,000 cm$^{-3}$)", ylabel="Altitude (km)",
               ylim=(150, 600), title="Electron density above the sounder")
        ax.grid(alpha=.25)
        ax.legend()
        profile_figure = (figure_output or output_prefix.with_name(
            output_prefix.name + "_density_profile.png"))
        profile_figure.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(profile_figure, dpi=250)
        plt.close(fig)
    summary_path = output_prefix.with_name(output_prefix.name + "_summary.json")
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    evaluation_path = output_prefix.with_name(output_prefix.name + "_evaluations.json")
    evaluation_path.write_text(json.dumps([
        {"candidate_index": index, "population": (index - 1) // 40,
         "density_scale": row["density_scale"],
         "hmf2_shift_km": row["hmf2_shift_km"],
         "f2_width_scale": row.get("f2_width_scale", 1.0),
         "score": row["score"]}
        for index, row in enumerate(evaluations, start=1)
    ], indent=2) + "\n")
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("state", type=Path)
    parser.add_argument("output_prefix", type=Path)
    parser.add_argument("--background", type=Path)
    parser.add_argument("--truth-density", type=Path)
    parser.add_argument("--figure-output", type=Path)
    args = parser.parse_args()
    summarize(args.state, args.output_prefix, args.background,
              args.truth_density, args.figure_output)


if __name__ == "__main__":
    main()
