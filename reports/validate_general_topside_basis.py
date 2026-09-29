"""Validate the IRI-2020 shape basis on held-out IRI-2016 vertical profiles.

These scalar O-mode ionograms check representational flexibility and inversion
conditioning. The separate saved 3-D O/X ray test checks transfer to PHaRLAP.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from python_raytrace.general_topside_inverse import (  # noqa: E402
    ChapmanF2, IRITopsideBasis, PLASMA_MHZ_PER_SQRT_CM3, fit_iri_basis_ionogram,
    fit_vertical_ionogram, _NODES, _WEIGHTS,
)

DATA = ROOT / "reports/data"
OUT = DATA / "general_topside_iri_case"
FREQUENCIES = np.round(np.arange(2.0, 10.0001, 0.1), 10)
SPACECRAFT_ALT_KM = 800.0


def tabulated_ionogram(altitude: np.ndarray, density: np.ndarray) -> np.ndarray:
    """Independent tabulated forward operator for a vertical scalar plasma."""
    peak = int(np.argmax(density))
    hm = float(altitude[peak])
    fof2 = PLASMA_MHZ_PER_SQRT_CM3 * np.sqrt(density[peak])
    interpolator = PchipInterpolator(altitude, np.log(density))
    rows = []
    for index, f in enumerate(FREQUENCIES):
        if f >= fof2:
            continue

        def q(z):
            return 1.0 - PLASMA_MHZ_PER_SQRT_CM3**2 * np.exp(interpolator(z)) / f**2

        if q(SPACECRAFT_ALT_KM) <= 0:
            continue
        turning = brentq(q, hm, SPACECRAFT_ALT_KM, xtol=1e-8)
        maximum = np.sqrt(SPACECRAFT_ALT_KM - turning)
        u = 0.5 * maximum * (_NODES + 1.0)
        group_range = maximum * np.dot(
            _WEIGHTS, 2.0 * u / np.sqrt(np.maximum(q(turning + u*u), 1e-14)))
        rows.append((index, 1, group_range, 0.0, 0.0))
    return np.asarray(rows)


def main() -> None:
    basis = IRITopsideBasis.read(DATA / "iri2020_global_topside_basis.npz")
    with np.load(DATA / "d_inverse_iri_truth_density.npz", allow_pickle=False) as source:
        lat = source["latitudes_deg"]
        lon = source["longitudes_deg"]
        altitude = source["altitudes_km"]
        density_grid = source["electron_density_cm3"]
    positions = [(i, j) for i in (6, 11, 16, 21) for j in (3, 10)]
    results = []
    for i, j in positions:
        truth = density_grid[i, j]
        records = tabulated_ionogram(altitude, truth)
        if len(records) < 15:
            raise ValueError(f"Too few returns at {lat[i]}, {lon[j]}")
        new = fit_iri_basis_ionogram(records, FREQUENCIES, SPACECRAFT_ALT_KM,
                                     basis, fit_x_mode=False)
        chapman = fit_vertical_ionogram(records, FREQUENCIES,
                                       SPACECRAFT_ALT_KM, fit_x_mode=False)
        hm = float(altitude[np.argmax(truth)])
        peak = float(np.max(truth))
        top = (altitude >= hm) & (altitude <= 600.0)

        def metrics(fit):
            predicted = fit.layer.density_cm3(altitude)
            return {"hmf2_error_km": float(fit.layer.hmf2_km - hm),
                    "fof2_error_mhz": float(fit.layer.fof2_mhz
                                            - PLASMA_MHZ_PER_SQRT_CM3 * np.sqrt(peak)),
                    "topside_normalized_rms_percent": float(
                        100 * np.sqrt(np.mean((predicted[top] - truth[top])**2)) / peak),
                    "range_mae_km": fit.range_mae_km[1]}

        results.append({"latitude_deg": float(lat[i]), "longitude_deg": float(lon[j]),
                        "truth_hmf2_km": hm, "truth_fof2_mhz": float(
                            PLASMA_MHZ_PER_SQRT_CM3 * np.sqrt(peak)),
                        "iri_basis": metrics(new), "chapman": metrics(chapman),
                        "selected_family_by_ridge_mae": (
                            "iri_basis" if new.range_mae_km[1] < chapman.range_mae_km[1]
                            else "chapman")})
    chapman_cases = []
    for layer in (ChapmanF2(5.4, 260.0, 55.0, 45.0),
                  ChapmanF2(5.4, 260.0, 55.0, 70.0),
                  ChapmanF2(5.4, 260.0, 55.0, 95.0),
                  ChapmanF2(6.5, 320.0, 55.0, 75.0, 0.5)):
        rows = []
        for index, frequency in enumerate(FREQUENCIES):
            group_range = layer.group_range_km(frequency, SPACECRAFT_ALT_KM)
            if np.isfinite(group_range):
                rows.append((index, 1, group_range, 0.0, 0.0))
        records = np.asarray(rows)
        shape_fit = fit_iri_basis_ionogram(records, FREQUENCIES,
                                           SPACECRAFT_ALT_KM, basis, fit_x_mode=False)
        chapman_fit = fit_vertical_ionogram(records, FREQUENCIES,
                                            SPACECRAFT_ALT_KM, fit_x_mode=False)
        truth_density = layer.density_cm3(altitude)
        truth_top = (altitude >= layer.hmf2_km) & (altitude <= 600.0)
        def analytic_metrics(fit):
            predicted = fit.layer.density_cm3(altitude)
            return {
                "fof2_error_mhz": float(fit.layer.fof2_mhz - layer.fof2_mhz),
                "hmf2_error_km": float(fit.layer.hmf2_km - layer.hmf2_km),
                "topside_normalized_rms_percent": float(
                    100 * np.sqrt(np.mean((predicted[truth_top] - truth_density[truth_top])**2))
                    / np.max(truth_density)),
                "range_mae_km": fit.range_mae_km[1],
            }
        chapman_cases.append({
            "truth_topside_scale_km": layer.topside_scale_km,
            "truth_tail_curvature": layer.tail_curvature,
            "iri_basis": analytic_metrics(shape_fit),
            "chapman": analytic_metrics(chapman_fit),
            "iri_basis_range_mae_km": shape_fit.range_mae_km[1],
            "chapman_range_mae_km": chapman_fit.range_mae_km[1],
            "selected_family_by_ridge_mae": (
                "iri_basis" if shape_fit.range_mae_km[1] < chapman_fit.range_mae_km[1]
                else "chapman"),
            "selected_hmf2_error_km": float(chapman_fit.layer.hmf2_km - layer.hmf2_km),
        })
    summary = {"basis_training_source": basis.source,
               "validation_density_source": "IRI-2016 1.11.1; eight separate positions",
               "forward_model": "vertical scalar-plasma group range on PCHIP of tabulated IRI-2016 density",
               "same_3d_ray_tracer_as_saved_ionograms": False,
               "limitation": "This test checks profile flexibility; it does not validate 3-D magnetoionic O/X rays.",
               "cases": results, "analytic_chapman_cases": chapman_cases}
    for model in ("iri_basis", "chapman"):
        summary[model + "_aggregate"] = {
            "hmf2_mae_km": float(np.mean([abs(row[model]["hmf2_error_km"]) for row in results])),
            "fof2_mae_mhz": float(np.mean([abs(row[model]["fof2_error_mhz"]) for row in results])),
            "topside_normalized_rms_mean_percent": float(np.mean([
                row[model]["topside_normalized_rms_percent"] for row in results])),
        }
    summary["family_selection"] = {
        "correct_iri2016_cases": sum(row["selected_family_by_ridge_mae"] == "iri_basis"
                                    for row in results),
        "correct_chapman_cases": sum(row["selected_family_by_ridge_mae"] == "chapman"
                                     for row in chapman_cases),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "scalar_validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    plot_validation(altitude, density_grid, positions, results, basis)
    print(json.dumps({k: v for k, v in summary.items()
                      if k.endswith("_aggregate") or k == "family_selection"}, indent=2))


def plot_validation(altitude: np.ndarray, density_grid: np.ndarray,
                    positions: list[tuple[int, int]], results: list[dict],
                    basis: IRITopsideBasis) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    i, j = positions[2]
    truth = density_grid[i, j]
    records = tabulated_ionogram(altitude, truth)
    iri_fit = fit_iri_basis_ionogram(records, FREQUENCIES, SPACECRAFT_ALT_KM,
                                     basis, fit_x_mode=False)
    chap_fit = fit_vertical_ionogram(records, FREQUENCIES,
                                    SPACECRAFT_ALT_KM, fit_x_mode=False)
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4.15), constrained_layout=True)
    axes[0].plot(truth / 1e5, altitude, color="black", linewidth=2.2,
                 label="IRI-2016 truth")
    axes[0].plot(iri_fit.layer.density_cm3(altitude) / 1e5, altitude,
                 color="#1261ac", linewidth=2, label="IRI basis fit")
    axes[0].plot(chap_fit.layer.density_cm3(altitude) / 1e5, altitude,
                 color="#a04a24", linewidth=2, linestyle="--",
                 label="Single Chapman fit")
    axes[0].set(xlabel="Electron density (100,000 cm$^{-3}$)",
                ylabel="Altitude (km)", ylim=(200, 600),
                title="Representative held-out profile")
    axes[0].grid(alpha=0.2)
    axes[0].legend(frameon=False)
    x = np.arange(len(results))
    axes[1].axhline(0, color="black", linewidth=0.9)
    axes[1].plot(x, [row["iri_basis"]["hmf2_error_km"] for row in results],
                 color="#1261ac", marker="o", label="IRI basis")
    axes[1].plot(x, [row["chapman"]["hmf2_error_km"] for row in results],
                 color="#a04a24", marker="s", linestyle="--", label="Single Chapman")
    axes[1].set(xlabel="Held-out IRI-2016 location", ylabel="hmF2 error (km)",
                xticks=x, xticklabels=[str(n+1) for n in x],
                title="Peak-height error across eight locations")
    axes[1].grid(alpha=0.2)
    axes[1].legend(frameon=False)
    fig.savefig(OUT / "scalar_validation.png", dpi=190)
    plt.close(fig)


if __name__ == "__main__":
    main()
