"""Check a density-anchored monotone-spline fit on the saved IRI wave pass.

This is a one-dimensional screening fit, not a replacement for 3-D full-ray
selection. ``prepare`` exposes only the synthetic in-situ density at 800 km;
``fit`` reads the saved O-mode ionograms, that scalar, and the previous
retrieval. Only ``evaluate`` opens the full independent truth profiles.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from python_raytrace.general_topside_inverse import (  # noqa: E402
    MonotoneSplineF2, PLASMA_MHZ_PER_SQRT_CM3, fit_spline_ionogram,
)

DATA = ROOT / "reports/data"
OUTPUT = DATA / "iri_wave_spline_local_density"
TRUTH = DATA / "lat_wave_pass_truth_density.npz"
BASELINE = DATA / "lat_wave_doppler_peak_round3/wide_density.npz"
OBSERVED = DATA / "lat_wave_final_truth_ionograms"
MANIFEST = DATA / "lat_wave_pass_manifest.json"


def positions() -> np.ndarray:
    rows = json.loads(MANIFEST.read_text())["profiles"]
    return np.array([(row["latitude_deg"], row["longitude_deg"]) for row in rows])


def sampled_density(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as source:
        altitude = np.asarray(source["altitudes_km"], dtype=float)
        interpolation = RegularGridInterpolator(
            (source["latitudes_deg"], source["longitudes_deg"]),
            source["electron_density_cm3"], bounds_error=True)
        profiles = np.asarray(interpolation(positions()), dtype=float)
    return altitude, profiles


def prepare() -> None:
    altitude, truth = sampled_density(TRUTH)
    measurements = np.array([np.interp(800.0, altitude, row) for row in truth])
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "local_density_observations.json").write_text(json.dumps({
        "source": "synthetic independent IRI-2016 wave density sampled at the 800 km sounder",
        "altitude_km": 800.0,
        "density_cm3_by_profile": measurements.tolist(),
        "full_truth_profile_exposed_to_fit": False,
    }, indent=2) + "\n")


def fit() -> None:
    observation = json.loads((OUTPUT / "local_density_observations.json").read_text())
    if observation["altitude_km"] != 800.0:
        raise ValueError("Expected the 800 km in-situ observation")
    measurements = np.asarray(observation["density_cm3_by_profile"], dtype=float)
    altitude, baseline = sampled_density(BASELINE)
    if len(measurements) != 20 or np.any(measurements <= 0):
        raise ValueError("Expected 20 positive local densities")
    fitted = np.empty_like(baseline)
    rows = []
    for index in range(1, 21):
        with np.load(OBSERVED / f"ionogram_{index:02d}.npz", allow_pickle=False) as ionogram:
            records = np.asarray(ionogram["records"], dtype=float)
            frequencies = np.asarray(ionogram["frequencies_mhz"], dtype=float)
            doppler = np.asarray(ionogram["spacecraft_doppler_hz"], dtype=float)
        # The largest ungated O nose can exceed the local foF2 when a ray
        # samples the wave crest away from the spacecraft. Retain the same
        # near-nadir Doppler gate used in the previous joint retrieval.
        records = records[(records[:, 1] == 1) & (abs(doppler) <= 15.0)]
        # The old ionogram-selected foF2 is an observational prior, not truth.
        anchor = PLASMA_MHZ_PER_SQRT_CM3 * np.sqrt(np.max(baseline[index - 1]))
        local_plasma = PLASMA_MHZ_PER_SQRT_CM3 * np.sqrt(measurements[index - 1])
        result = fit_spline_ionogram(
            records, frequencies, 800.0, regularize_tail=True,
            fof2_anchor_mhz=float(anchor),
            spacecraft_plasma_mhz_anchor=float(local_plasma))
        layer: MonotoneSplineF2 = result.layer
        fitted[index - 1] = layer.density_cm3(altitude)
        rows.append({
            "index": index, "fof2_mhz": layer.fof2_mhz,
            "hmf2_km": layer.hmf2_km,
            "bottomside_scale_km": layer.bottomside_scale_km,
            "log_slope_per_km": layer.log_slope_per_km.tolist(),
            "baseline_fof2_anchor_mhz": float(anchor),
            "local_density_observation_cm3": float(measurements[index - 1]),
            "fitted_800km_density_cm3": float(layer.density_cm3(800.0)),
            "one_dimensional_o_ridge_mae_km": result.range_mae_km[1],
        })
        print(f"fitted profile {index:02d}: foF2={layer.fof2_mhz:.3f} MHz, "
              f"hmF2={layer.hmf2_km:.1f} km", flush=True)
    np.savez_compressed(OUTPUT / "fitted_profiles.npz", altitude_km=altitude,
                        latitudes_deg=positions()[:, 0],
                        electron_density_cm3=fitted)
    (OUTPUT / "fits.json").write_text(json.dumps({
        "truth_profile_read_during_fit": False,
        "full_3d_ray_ionogram_selection_complete": False,
        "method": "independent eight-slope monotone topside fits to 15 Hz gated O ridges; "
                  "prior full-ray foF2 and synthetic local Ne at 800 km anchored",
        "profiles": rows,
    }, indent=2) + "\n")


def tail() -> None:
    """Anchor the previous retrieval at 800 km without changing its F2 peak."""
    observation = json.loads((OUTPUT / "local_density_observations.json").read_text())
    measured = np.asarray(observation["density_cm3_by_profile"], dtype=float)
    altitude, baseline = sampled_density(BASELINE)
    local_baseline = np.array([np.interp(800.0, altitude, row) for row in baseline])
    if np.any(local_baseline <= 0) or np.any(measured <= 0):
        raise ValueError("Nonpositive local density")
    # Smooth in altitude, with no modification through 600 km. This is a
    # candidate for full-ray checking, not an ionogram-selected result.
    x = np.clip((altitude - 600.0) / 200.0, 0.0, 1.0)
    weight = x * x * (3.0 - 2.0 * x)
    ratio = measured / local_baseline
    revised = baseline * np.exp(np.log(ratio[:, None]) * weight[None, :])
    if np.any(np.diff(revised[:, altitude >= 300], axis=1) > 0):
        raise ValueError("Tail adjustment destroyed monotone topside")
    np.savez_compressed(OUTPUT / "tail_anchored_profiles.npz",
                        altitude_km=altitude, latitudes_deg=positions()[:, 0],
                        electron_density_cm3=revised)
    (OUTPUT / "tail_candidate.json").write_text(json.dumps({
        "truth_profile_read_during_fit": False,
        "full_3d_ray_ionogram_selection_complete": False,
        "method": "smooth log-density taper from 600 to 800 km, exact local-density anchor above 800 km",
        "density_ratio_at_800km_by_profile": ratio.tolist(),
        "peak_profile_unchanged": True,
    }, indent=2) + "\n")


def evaluate() -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fit_record = json.loads((OUTPUT / "fits.json").read_text())
    if fit_record["truth_profile_read_during_fit"]:
        raise ValueError("Fit opened the truth density")
    altitude, truth = sampled_density(TRUTH)
    base_altitude, baseline = sampled_density(BASELINE)
    with np.load(OUTPUT / "fitted_profiles.npz", allow_pickle=False) as source:
        fit_altitude = source["altitude_km"]
        fitted = source["electron_density_cm3"]
    with np.load(OUTPUT / "tail_anchored_profiles.npz", allow_pickle=False) as source:
        tail_altitude = source["altitude_km"]
        tail_fit = source["electron_density_cm3"]
    if not (np.array_equal(altitude, base_altitude)
            and np.array_equal(altitude, fit_altitude)
            and np.array_equal(altitude, tail_altitude)):
        raise ValueError("Altitude axes differ")
    mask = (altitude >= 150) & (altitude <= 600)
    upper_mask = (altitude >= 300) & (altitude <= 800)
    true_peak = truth.max(axis=1)
    true_height = altitude[np.argmax(truth, axis=1)]

    def metrics(profiles: np.ndarray) -> dict:
        peak = profiles.max(axis=1)
        height = altitude[np.argmax(profiles, axis=1)]
        return {
            "peak_density_mae_percent": float(np.mean(abs(100 * (peak / true_peak - 1)))),
            "fof2_mae_mhz": float(np.mean(abs(PLASMA_MHZ_PER_SQRT_CM3
                                                  * (np.sqrt(peak) - np.sqrt(true_peak))))),
            "hmf2_mae_km": float(np.mean(abs(height - true_height))),
            "density_nrmse_150_to_600": float(np.mean(
                np.sqrt(np.mean((profiles[:, mask] - truth[:, mask]) ** 2, axis=1))
                / true_peak)),
            "density_nrmse_300_to_800": float(np.mean(
                np.sqrt(np.mean((profiles[:, upper_mask] - truth[:, upper_mask]) ** 2,
                                axis=1)) / true_peak)),
            "local_800km_density_mae_percent": float(np.mean(abs(100 * (
                profiles[:, np.argmin(abs(altitude - 800))]
                / truth[:, np.argmin(abs(altitude - 800))] - 1)))),
            "hmf2_error_km_by_profile": (height - true_height).tolist(),
            "peak_density_error_percent_by_profile": (100 * (peak / true_peak - 1)).tolist(),
        }

    result = {
        "case": "20-position independent IRI-2016 latitude wave",
        "candidate_status": "one-dimensional O-mode screening only; no 3-D full-ray selection",
        "synthetic_local_density_observation": True,
        "baseline": metrics(baseline),
        "spline": metrics(fitted),
        "tail_anchor": metrics(tail_fit),
        "spline_positions_missing_local_density_by_over_5_percent": [
            row["index"] for row in fit_record["profiles"]
            if abs(row["fitted_800km_density_cm3"]
                   / row["local_density_observation_cm3"] - 1) > .05],
    }
    (OUTPUT / "evaluation.json").write_text(json.dumps(result, indent=2) + "\n")
    latitude = positions()[:, 0]
    figure, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for row, index in enumerate((8, 14)):
        ax = axes[row, 0]
        ax.plot(truth[index - 1] / 1e5, altitude, "k", label="IRI wave truth")
        ax.plot(baseline[index - 1] / 1e5, altitude, color="#777777",
                label="Previous full-ray retrieval")
        ax.plot(fitted[index - 1] / 1e5, altitude, color="#257e56",
                label="Spline screen (soft 800 km anchor)")
        ax.plot(tail_fit[index - 1] / 1e5, altitude, color="#c06b25",
                linestyle=":", label="Conservative 800 km anchor")
        ax.set(xlim=(0, 5), ylim=(150, 850),
               xlabel="Density ($10^5$ cm$^{-3}$)", ylabel="Altitude (km)",
               title=f"Profile {index:02d}")
        ax.grid(alpha=.2)
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].plot(latitude, result["baseline"]["hmf2_error_km_by_profile"],
                    "o--", color="#777777", label="Previous")
    axes[0, 1].plot(latitude, result["spline"]["hmf2_error_km_by_profile"],
                    "o-", color="#257e56", label="Spline screen")
    axes[0, 1].set(ylabel="hmF2 error (km)", title="Peak height")
    axes[1, 1].plot(latitude, result["baseline"]["peak_density_error_percent_by_profile"],
                    "o--", color="#777777")
    axes[1, 1].plot(latitude, result["spline"]["peak_density_error_percent_by_profile"],
                    "o-", color="#257e56")
    axes[1, 1].set(xlabel="Latitude (degrees)", ylabel="Peak-density error (%)",
                   title="Peak density")
    for ax in axes[:, 1]:
        ax.axhline(0, color="black", linewidth=.8)
        ax.grid(alpha=.2)
    figure.savefig(OUTPUT / "spline_screen.png", dpi=180)
    plt.close(figure)
    print(json.dumps({key: {m: v for m, v in value.items() if "by_profile" not in m}
                      if isinstance(value, dict) else value
                      for key, value in result.items()}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("prepare", "fit", "tail", "evaluate"))
    arguments = parser.parse_args()
    {"prepare": prepare, "fit": fit, "tail": tail, "evaluate": evaluate}[arguments.stage]()
