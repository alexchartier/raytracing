"""Build density cuts and peak diagnostics from the frozen vertical and oblique fits."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "reports/data"
FIGURES = ROOT / "reports/figures"
TRUTH = DATA / "lat_wave_pass_truth_density.npz"
VERTICAL = DATA / "lat_wave_doppler_peak_round3/wide_density.npz"
OBLIQUE = DATA / "oblique_wave_600km/gain0_sigma100_height1_spline/density.npz"
PEAK_EVALUATION = DATA / "vertical_oblique_peak_evaluation.json"


def density_section(path: Path, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        altitude = np.asarray(data["altitudes_km"], dtype=float)
        interpolate = RegularGridInterpolator(
            (data["latitudes_deg"], data["longitudes_deg"]),
            data["electron_density_cm3"], bounds_error=True)
        density = np.asarray(interpolate(points), dtype=float)
    return altitude, density


def altitude_cut(path: Path, latitude: float, longitude: float) -> tuple[np.ndarray, np.ndarray]:
    altitude, density = density_section(path, np.array([[latitude, longitude]]))
    return altitude, density[0]


def f2_peaks(altitude: np.ndarray, density: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Estimate foF2 and hmF2 with a 3-point quadratic about each sampled maximum."""
    if density.ndim != 2 or density.shape[1] != len(altitude):
        raise ValueError("Expected one density profile per row")
    spacing = np.diff(altitude)
    if not np.allclose(spacing, 20.0):
        raise ValueError("Expected the saved 20 km altitude grid")
    indices = np.argmax(density, axis=1)
    if np.any((indices == 0) | (indices == len(altitude) - 1)):
        raise ValueError("F2 peak lies on an altitude-grid boundary")
    rows = np.arange(len(indices))
    before = density[rows, indices - 1]
    center = density[rows, indices]
    after = density[rows, indices + 1]
    curvature = before - 2 * center + after
    if np.any(curvature >= 0):
        raise ValueError("Expected concave local peaks")
    shift = 0.5 * (before - after) / curvature * spacing[0]
    if np.any(np.abs(shift) > spacing[0]):
        raise ValueError("Quadratic vertex lies beyond its three samples")
    peak_density = center - (after - before) ** 2 / (8 * curvature)
    return 0.00898 * np.sqrt(peak_density), altitude[indices] + shift


def peak_metrics(truth: np.ndarray, retrieved: np.ndarray) -> dict[str, float]:
    residual = retrieved - truth
    return {"mae": float(np.mean(np.abs(residual))),
            "bias": float(np.mean(residual)),
            "max_absolute_error": float(np.max(np.abs(residual)))}


def evaluate_peaks(vertical_points: np.ndarray) -> None:
    midpoint_points = []
    for index in range(1, 21):
        with np.load(DATA / "oblique_wave_600km/truth" /
                     f"ionogram_{index:02d}.npz", allow_pickle=False) as link:
            midpoint_points.append([
                (float(link["tx_lat_deg"]) + float(link["rx_lat_deg"])) / 2,
                float(link["tx_lon_deg"]),
            ])
    midpoint_points = np.asarray(midpoint_points)
    cases = (
        ("vertical", vertical_points, VERTICAL),
        ("oblique_start", midpoint_points, VERTICAL),
        ("oblique_selected", midpoint_points, OBLIQUE),
    )
    result = {
        "method": "Three-point quadratic around each maximum after horizontal interpolation; foF2 = 0.00898 sqrt(N_peak in cm^-3) MHz",
        "altitude_grid_spacing_km": 20.0,
        "truth_density": str(TRUTH.relative_to(ROOT)),
        "evaluations": {},
    }
    curves = {}
    for name, points, retrieval in cases:
        altitude, truth_density = density_section(TRUTH, points)
        retrieval_altitude, retrieved_density = density_section(retrieval, points)
        if not np.array_equal(altitude, retrieval_altitude):
            raise ValueError("Truth and retrieval altitude grids differ")
        true_fof2, true_hmf2 = f2_peaks(altitude, truth_density)
        fit_fof2, fit_hmf2 = f2_peaks(altitude, retrieved_density)
        curves[name] = (points[:, 0], true_fof2, fit_fof2, true_hmf2, fit_hmf2)
        result["evaluations"][name] = {
            "locations": ("sounder_positions" if name == "vertical"
                          else "oblique_link_midpoints"),
            "retrieved_density": str(retrieval.relative_to(ROOT)),
            "fof2_mhz": peak_metrics(true_fof2, fit_fof2),
            "hmf2_km": peak_metrics(true_hmf2, fit_hmf2),
            "latitude_deg": points[:, 0].tolist(),
            "truth_fof2_mhz": true_fof2.tolist(),
            "retrieved_fof2_mhz": fit_fof2.tolist(),
            "truth_hmf2_km": true_hmf2.tolist(),
            "retrieved_hmf2_km": fit_hmf2.tolist(),
        }
    PEAK_EVALUATION.write_text(json.dumps(result, indent=2) + "\n")

    fig, axes = plt.subplots(2, 2, figsize=(11, 6.5), sharey="row",
                             constrained_layout=True)
    for column, (name, label) in enumerate((("vertical", "Vertical sounders"),
                                            ("oblique_selected", "Oblique link midpoints"))):
        lat, true_fof2, fit_fof2, true_hmf2, fit_hmf2 = curves[name]
        for row, (truth, retrieved, unit, key) in enumerate((
                (true_fof2, fit_fof2, "MHz", "fof2_mhz"),
                (true_hmf2, fit_hmf2, "km", "hmf2_km"))):
            ax = axes[row, column]
            ax.plot(lat, truth, "o-", color="#242424", ms=3.7,
                    lw=1.6, label="Independent truth")
            ax.plot(lat, retrieved, "s-", color="#208054", ms=3.7,
                    lw=1.6, label="Selected retrieval")
            mae = result["evaluations"][name][key]["mae"]
            ax.set_title(f"{label}  |  MAE {mae:.3f} {unit}" if row == 0 else
                         f"{label}  |  MAE {mae:.1f} {unit}", fontsize=11)
            ax.grid(alpha=0.2)
            ax.set_xlabel("Latitude (degrees)")
        axes[0, column].legend(fontsize=8, loc="upper left")
    axes[0, 0].set_ylabel("foF2 (MHz)")
    axes[1, 0].set_ylabel("hmF2 (km)")
    fig.savefig(FIGURES / "vertical_oblique_peak_frequency_height.png", dpi=220)
    plt.close(fig)


def plot_cut(latitude: float, longitude: float, retrieval: Path,
             title: str, output: Path) -> None:
    altitude, truth = altitude_cut(TRUTH, latitude, longitude)
    retrieved_altitude, retrieved = altitude_cut(retrieval, latitude, longitude)
    if not np.array_equal(altitude, retrieved_altitude):
        raise ValueError("Truth and retrieval altitude grids differ")

    mask = (altitude >= 150) & (altitude <= 600)
    true_peak = np.max(truth)
    fit_peak = np.max(retrieved)
    peak_error = 100 * (fit_peak / true_peak - 1)
    profile_error = 100 * np.sqrt(np.mean((retrieved[mask] - truth[mask]) ** 2)) / true_peak

    fig, ax = plt.subplots(figsize=(7.4, 3.2), constrained_layout=True)
    ax.plot(truth[mask] / 1e5, altitude[mask], color="#242424", lw=2.2,
            label="Independent synthetic truth")
    ax.plot(retrieved[mask] / 1e5, altitude[mask], color="#208054", lw=2.2,
            label="Selected retrieval")
    ax.set(xlabel=r"Electron density ($10^5$ cm$^{-3}$)",
           ylabel="Altitude (km)", ylim=(150, 600), title=title)
    ax.grid(alpha=0.2)
    ax.legend(loc="upper right", fontsize=9)
    ax.text(0.98, 0.06,
            f"Peak error {peak_error:+.2f}%   |   Profile NRMSE {profile_error:.2f}%",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "#dddddd", "alpha": 0.9})
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    vertical_choice = json.loads((DATA / "lat_wave_doppler_peak_round3/final_selection.json").read_text())
    oblique_choice = json.loads((DATA / "oblique_wave_600km/selection.json").read_text())
    if (vertical_choice["chosen_candidate"] != "wide"
            or oblique_choice["chosen_candidate"] != "gain0_sigma100_height1_spline"):
        raise ValueError("The frozen selections have changed; update the report inputs")

    manifest = json.loads((DATA / "lat_wave_pass_manifest.json").read_text())
    row = manifest["profiles"][13]
    FIGURES.mkdir(exist_ok=True)
    plot_cut(float(row["latitude_deg"]), float(row["longitude_deg"]), VERTICAL,
             f"Vertical profile 14 ({float(row['latitude_deg']):.1f}°)",
             FIGURES / "vertical_oblique_vertical_profile14.png")

    with np.load(DATA / "oblique_wave_600km/truth/ionogram_14.npz", allow_pickle=False) as link:
        latitude = (float(link["tx_lat_deg"]) + float(link["rx_lat_deg"])) / 2
        longitude = float(link["tx_lon_deg"])
    plot_cut(latitude, longitude, OBLIQUE,
             f"Oblique link 14 midpoint ({latitude:.1f}°)",
             FIGURES / "vertical_oblique_oblique_link14.png")
    vertical_points = np.array([[row["latitude_deg"], row["longitude_deg"]]
                                for row in manifest["profiles"]], dtype=float)
    evaluate_peaks(vertical_points)


if __name__ == "__main__":
    main()
