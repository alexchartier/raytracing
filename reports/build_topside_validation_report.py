"""Reproduce a small, clearly labeled topside retrieval diagnostic PDF.

Run from the repository root with ``python3 reports/build_topside_validation_report.py``.
Requires the installed PyIRI/PHaRLAP runtime. No AMPERE data file is required.
"""

from __future__ import annotations

import datetime as dt
import json
import math
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from python_raytrace.multisat_topside_inverse_demo import (  # noqa: E402
    IonosphereFitParams,
    SyntheticDataset,
    TopsideInverseConfig,
    _apply_fit_params_to_grid,
    _gmst_from_jd,
    _render_case_observables,
    _subset_problem_cases,
    build_inverse_problem,
    dataset_cost,
    simulate_dataset,
)


def write_orbit_fixture(path: Path) -> None:
    """Generate reproducible Earth orbit geometry, not an AMPERE observation."""
    hours = np.arange(11.5, 12.5001, 1.0 / 60.0)
    radius_km = 7178.0
    inclination_rad = math.radians(55.0)
    positions_m = []
    for hour in hours:
        when = dt.datetime(2010, 1, 1, tzinfo=dt.timezone.utc) + dt.timedelta(hours=float(hour))
        jd = when.timestamp() / 86400.0 + 2440587.5
        phase = (float(hour) - 12.0) * 2.0 * math.pi / 1.6 + _gmst_from_jd(jd)
        positions_m.append(
            radius_km * 1000.0 * np.array(
                [math.cos(phase), math.sin(phase) * math.cos(inclination_rad), math.sin(phase) * math.sin(inclination_rad)]
            )
        )
    with netCDF4.Dataset(path, "w") as dataset:
        dataset.createDimension("obs", len(hours))
        dataset.createDimension("xyz", 3)
        for name, values, dims, dtype in (
            ("time", hours, ("obs",), "f8"),
            ("plane_num", np.ones(len(hours), dtype=int), ("obs",), "i4"),
            ("pseudo_sv_num", np.full(len(hours), 121, dtype=int), ("obs",), "i4"),
            ("pos_eci", np.asarray(positions_m), ("obs", "xyz"), "f8"),
            ("year", np.full(len(hours), 2010, dtype=int), ("obs",), "i4"),
            ("doy", np.ones(len(hours), dtype=int), ("obs",), "i4"),
        ):
            dataset.createVariable(name, dtype, dims)[:] = values


def dataset_from_case(problem, case_observables) -> SyntheticDataset:
    return SyntheticDataset(
        frequencies_mhz=problem.frequencies_mhz,
        range_edges_km=problem.range_edges_km,
        range_centers_km=problem.range_centers_km,
        cases=(case_observables,),
    )


def main() -> None:
    output = ROOT / "reports" / "topside_sounder_status_2026-09-25.pdf"
    metrics_path = ROOT / "reports" / "topside_validation_metrics.json"
    true_scale = 1.12
    scales = np.array([0.80, 0.90, 1.00, 1.10, 1.20, 1.30], dtype=float)
    zero_wave = IonosphereFitParams(
        density_scale=1.0,
        hmf2_shift_km=0.0,
        wave_amplitude_fraction=0.0,
        wave_phase_rad=0.0,
        wave_bearing_deg=0.0,
    )

    with tempfile.TemporaryDirectory(prefix="topside_validation_") as temporary_dir:
        temporary = Path(temporary_dir)
        orbit_path = temporary / "generated_orbit.nc"
        write_orbit_fixture(orbit_path)
        config = replace(
            TopsideInverseConfig(),
            ampere_file=orbit_path,
            planes=(1,),
            frequencies_mhz=(4.0, 5.0, 6.0),
            vertical_elevation_count=8,
            vertical_azimuth_step_deg=90.0,
            oblique_elevation_count=9,
            oblique_bearing_count=7,
            grid_lat_step_deg=2.0,
            grid_lon_step_deg=4.0,
            grid_alt_step_km=20.0,
            d_region_model="none",
            grid_cache_path=temporary / "background.nc",
        )
        problem = _subset_problem_cases(build_inverse_problem(config), ["plane01_sv121_vertical"])
        case = problem.cases[0]
        background = problem.background_grids[0]

        # Closed-model experiment: the retrieval family also generates truth.
        same_model_truth = simulate_dataset(problem, replace(zero_wave, density_scale=true_scale))

        # Structural mismatch experiment: a separate altitude-dependent formula
        # generates density. The ray tracer and PyIRI background remain shared.
        heights = np.asarray(background.altitudes_km, dtype=float)
        altitude_factor = 1.10 + 0.10 * np.exp(-0.5 * ((heights - 320.0) / 55.0) ** 2)
        independent_grid = replace(
            background,
            iono_en_grid=np.asarray(background.iono_en_grid) * altitude_factor[None, None, :],
            iono_en_grid_5=np.asarray(background.iono_en_grid_5) * altitude_factor[None, None, :],
        )
        structural_truth = dataset_from_case(
            problem, _render_case_observables(problem, case, independent_grid, None)
        )

        no_return_factor = 0.92 + 0.33 * np.exp(-0.5 * ((heights - 320.0) / 55.0) ** 2)
        no_return_grid = replace(
            background,
            iono_en_grid=np.asarray(background.iono_en_grid) * no_return_factor[None, None, :],
            iono_en_grid_5=np.asarray(background.iono_en_grid_5) * no_return_factor[None, None, :],
        )
        no_return_truth = dataset_from_case(
            problem, _render_case_observables(problem, case, no_return_grid, None)
        )
        if float(np.sum(same_model_truth.cases[0].total_image)) <= 0.0:
            raise RuntimeError("closed-model truth has no detectable ionogram returns")
        if float(np.sum(structural_truth.cases[0].total_image)) <= 0.0:
            raise RuntimeError("structural-mismatch truth has no detectable ionogram returns")

        predictions = [simulate_dataset(problem, replace(zero_wave, density_scale=float(scale))) for scale in scales]
        closed_costs = np.array([dataset_cost(same_model_truth, prediction) for prediction in predictions])
        structural_costs = np.array([dataset_cost(structural_truth, prediction) for prediction in predictions])
        closed_index = int(np.argmin(closed_costs))
        structural_index = int(np.argmin(structural_costs))

        fit_grid = _apply_fit_params_to_grid(problem, background, replace(zero_wave, density_scale=float(scales[structural_index])))
        lat_index = len(background.latitudes_deg) // 2
        lon_index = len(background.longitudes_deg) // 2
        truth_profile = independent_grid.iono_en_grid[lat_index, lon_index, :]
        fit_profile = fit_grid.iono_en_grid[lat_index, lon_index, :]
        profile_mask = (heights >= 150.0) & (heights <= 450.0)
        profile_nrmse = float(
            np.sqrt(np.mean((truth_profile[profile_mask] - fit_profile[profile_mask]) ** 2))
            / np.sqrt(np.mean(truth_profile[profile_mask] ** 2))
        )
        no_return_profile = no_return_grid.iono_en_grid[lat_index, lon_index, :]
        zero_fit_profile = _apply_fit_params_to_grid(
            problem, background, replace(zero_wave, density_scale=float(scales[0]))
        ).iono_en_grid[lat_index, lon_index, :]
        no_return_nrmse = float(
            np.sqrt(np.mean((no_return_profile[profile_mask] - zero_fit_profile[profile_mask]) ** 2))
            / np.sqrt(np.mean(no_return_profile[profile_mask] ** 2))
        )

        metrics = {
            "scenario": "generated Earth orbit; PyIRI background; one vertical case; 4, 5, 6 MHz; density-scale retrieval only",
            "same_model": {
                "true_scale": true_scale,
                "retrieved_scale": float(scales[closed_index]),
                "absolute_scale_error": float(abs(scales[closed_index] - true_scale)),
                "objective": float(closed_costs[closed_index]),
                "generated_ionogram_power_sum": float(np.sum(same_model_truth.cases[0].total_image)),
            },
            "structural_mismatch": {
                "truth_formula": "PyIRI density times [1.10 + 0.10 exp(-0.5 ((altitude_km - 320)/55)^2)]",
                "retrieved_scale": float(scales[structural_index]),
                "objective": float(structural_costs[structural_index]),
                "profile_nrmse_150_450_km": profile_nrmse,
                "generated_ionogram_power_sum": float(np.sum(structural_truth.cases[0].total_image)),
                "independent_background": False,
                "independent_ray_tracer": False,
            },
            "no_return_control": {
                "generated_ionogram_power_sum": float(np.sum(no_return_truth.cases[0].total_image)),
                "zero_return_candidate_scale": float(scales[0]),
                "objective": float(dataset_cost(no_return_truth, predictions[0])),
                "profile_nrmse_150_450_km": no_return_nrmse,
            },
        }
        metrics_path.write_text(json.dumps(metrics, indent=2) + "\n")

        with PdfPages(output, metadata={"Title": "Topside sounder retrieval status and validation audit"}) as pdf:
            fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.5), constrained_layout=True)
            ax = axes[0, 0]
            ax.plot(scales, closed_costs, "o-", color="#176B87")
            ax.axvline(true_scale, color="#B6483C", linestyle="--", label=f"Generating scale = {true_scale:.2f}")
            ax.axvline(scales[closed_index], color="#1A774E", linestyle=":", label=f"Retrieved = {scales[closed_index]:.2f}")
            ax.set(title="A  Closed-model, coarse grid search", xlabel="Density scale", ylabel="Ionogram objective")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.25)

            ax = axes[0, 1]
            ax.plot(scales, structural_costs, "o-", color="#7553A5")
            ax.axvline(scales[structural_index], color="#1A774E", linestyle=":", label=f"Retrieved = {scales[structural_index]:.2f}")
            ax.set(title="B  Separate profile-shape formula", xlabel="Fitted density scale", ylabel="Ionogram objective")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.25)

            ax = axes[1, 0]
            ax.plot(truth_profile[profile_mask], heights[profile_mask], label="Generated profile", color="#B6483C")
            ax.plot(fit_profile[profile_mask], heights[profile_mask], label="Retrieved profile", color="#176B87")
            ax.set(title=f"C  Profile NRMSE = {profile_nrmse:.1%} (150–450 km)", xlabel="Electron density (cm$^{-3}$)", ylabel="Altitude (km)")
            ax.legend(fontsize=8)
            ax.grid(alpha=0.25)

            ax = axes[1, 1]
            observed_image = structural_truth.cases[0].total_image
            fitted_image = predictions[structural_index].cases[0].total_image
            image_difference = fitted_image - observed_image
            vmax = max(float(np.max(np.abs(image_difference))), 1e-6)
            image = ax.imshow(
                image_difference.T,
                extent=(3.5, 6.5, float(problem.range_edges_km[-1]), float(problem.range_edges_km[0])),
                aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax,
            )
            ax.set(title="D  Fitted minus generated ionogram", xlabel="Frequency (MHz)", ylabel="Virtual range (km)")
            fig.colorbar(image, ax=ax, label="Normalized power difference")
            fig.suptitle("Topside retrieval diagnostic — generated Earth case, not Mars validation", fontsize=15)
            pdf.savefig(fig)
            plt.close(fig)

            fig = plt.figure(figsize=(8.5, 11))
            ax = fig.add_axes([0.09, 0.07, 0.82, 0.86]); ax.axis("off")
            lines = [
                ("Topside retrieval: evidence and status", 16, True),
                ("25 September 2026 | raytracing / python_raytrace", 9, False),
                ("What the figures establish", 12, True),
                (f"Closed-model test: true density scale {true_scale:.2f}; six-point grid search returned {scales[closed_index]:.2f}; absolute scale error {abs(scales[closed_index]-true_scale):.2f}; objective {closed_costs[closed_index]:.5g}. The true value is off-grid.", 10, False),
                (f"Separate profile-shape test: best scale {scales[structural_index]:.2f}; 150–450 km profile NRMSE {profile_nrmse:.1%}; ionogram objective {structural_costs[structural_index]:.5g}. This tests model mismatch, not a Mars measurement.", 10, False),
                (f"No-return control: both truth and a scale-{scales[0]:.2f} candidate have zero ionogram power, giving objective zero despite {no_return_nrmse:.1%} profile NRMSE. A zero objective alone is therefore not an accuracy claim.", 10, False),
                ("Truth independence audit", 12, True),
                ("The existing demo generates truth with simulate_dataset and fits it with simulate_dataset using the same PyIRI background, parameterized density field, and PHaRLAP ray tracer. Its truth is therefore not independent of retrieval assumptions.", 10, False),
                ("The new structural check generates the density profile with a separate altitude formula, outside the fitted scalar family. It still shares PyIRI, ray tracing, image formation, and generated Earth orbit geometry. It is only partially independent.", 10, False),
                ("No independent empirical Mars truth or measured Mars ionogram has yet been run through the retrieval. Accordingly, this report does not claim real-world retrieval accuracy.", 10, False),
                ("Public Mars reference candidate", 12, True),
                ("Huang et al. (2021) publish Mathematica expressions for 12 dayside ion species over 150–450 km, driven by EUV index, solar zenith angle, and magnetic elevation. Charge-balanced total ion density is a possible external electron-density reference, pending unit and license checks.", 10, False),
                ("github.com/hamsterping/Empirical-models-of-ion-density-distribution-in-the-dayside-Martian-ionosphere", 8.5, False),
                ("doi.org/10.1029/2021JA029226", 8.5, False),
                ("NeMars is a closer MARSIS electron-density model (doi.org/10.1016/j.icarus.2013.03.021), but I did not locate runnable public code. No public Mars-IRI or Mars-NeQuick implementation was verified.", 10, False),
                ("Next validation gate", 12, True),
                ("Drive observations from an external Mars density grid or measured MARSIS profiles, with Mars geometry and magnetic/neutral assumptions, then compare retrieved density to withheld truth on a common altitude and location grid. Report bias, NRMSE, and uncertainty across multiple cases.", 10, False),
            ]
            y = 0.98
            import textwrap
            for content, size, bold in lines:
                width = 96 if size >= 10 else 105
                chunks = textwrap.wrap(content, width=width, break_long_words=False, break_on_hyphens=False)
                if not chunks:
                    chunks = [""]
                for chunk in chunks:
                    ax.text(0, y, chunk, transform=ax.transAxes, va="top", fontsize=size,
                            fontweight="bold" if bold else "normal", color="#173449")
                    y -= 0.028 if size >= 12 else 0.021
                y -= 0.012 if bold else 0.015
            pdf.savefig(fig)
            plt.close(fig)

    print(json.dumps(metrics, indent=2))
    print(output)


if __name__ == "__main__":
    main()
