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
import textwrap
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

        # Render a denser frequency grid solely for the report. The selected
        # parameters above still come from the three-frequency search.
        display_frequencies = np.arange(4.0, 6.0001, 0.25)
        display_problem = replace(problem, frequencies_mhz=display_frequencies)
        closed_truth_display = simulate_dataset(display_problem, replace(zero_wave, density_scale=true_scale))
        closed_fit_display = simulate_dataset(
            display_problem, replace(zero_wave, density_scale=float(scales[closed_index]))
        )
        shape_truth_display = dataset_from_case(
            display_problem, _render_case_observables(display_problem, case, independent_grid, None)
        )
        shape_fit_display = simulate_dataset(
            display_problem, replace(zero_wave, density_scale=float(scales[structural_index]))
        )
        display_fields = tuple(
            data.cases[0].total_image
            for data in (closed_truth_display, closed_fit_display, shape_truth_display, shape_fit_display)
        )
        visible_range_indices = np.flatnonzero(np.max(np.stack(display_fields), axis=(0, 1)) > 0.01)
        if visible_range_indices.size:
            display_range_min = max(float(problem.range_edges_km[0]),
                                    float(problem.range_edges_km[visible_range_indices[0]]) - 100.0)
            display_range_max = min(float(problem.range_edges_km[-1]),
                                    float(problem.range_edges_km[visible_range_indices[-1] + 1]) + 100.0)
        else:
            display_range_min = float(problem.range_edges_km[0])
            display_range_max = float(problem.range_edges_km[-1])

        def ionogram_pair(figure, truth_field, fit_field, *, bottom=0.43, height=0.30):
            axes = (figure.add_axes([0.11, bottom, 0.32, height]),
                    figure.add_axes([0.55, bottom, 0.32, height]))
            image = None
            for ax, field, title in zip(axes, (truth_field, fit_field),
                                        ("Synthetic truth ionogram", "Retrieved ionogram")):
                image = ax.imshow(
                    field.T, origin="lower", aspect="auto", cmap="magma", vmin=0.0, vmax=1.0,
                    extent=(3.875, 6.125, float(problem.range_edges_km[0]), float(problem.range_edges_km[-1])),
                    interpolation="nearest",
                )
                ax.set_ylim(display_range_max, display_range_min)
                ax.set_xlim(3.875, 6.125)
                ax.set_title(title, fontsize=10, weight="bold")
                ax.set_xlabel("Frequency (MHz)")
                ax.set_ylabel("Virtual range (km)")
            colorbar_ax = figure.add_axes([0.89, bottom, 0.015, height])
            figure.colorbar(image, cax=colorbar_ax, label="Normalized power").ax.tick_params(labelsize=8)

        def page(title: str, number: int):
            figure = plt.figure(figsize=(8.5, 11), facecolor="white")
            figure.text(0.09, 0.948, title, fontsize=18, weight="bold", color="#173449", va="top")
            figure.text(0.09, 0.914, "TOPSIDE SOUNDER RETRIEVAL  |  25 SEPTEMBER 2026", fontsize=8.5,
                        color="#5B7180", va="top")
            figure.add_artist(plt.Line2D([0.09, 0.91], [0.895, 0.895], transform=figure.transFigure,
                                         color="#A8BBC6", linewidth=0.8))
            figure.text(0.09, 0.043, "raytracing  ·  generated Earth diagnostic  ·  no external truth data",
                        fontsize=8, color="#5B7180")
            figure.text(0.91, 0.043, f"{number} / 6", ha="right", fontsize=8, color="#5B7180")
            return figure

        def paragraph(figure, text: str, y: float, *, x: float = 0.09, width: int = 93,
                      size: float = 10.5, line_height: float = 0.023, color: str = "#173449") -> float:
            for row in textwrap.wrap(text, width=width, break_long_words=False, break_on_hyphens=False):
                figure.text(x, y, row, va="top", fontsize=size, color=color)
                y -= line_height
            return y

        def heading(figure, text: str, y: float) -> float:
            figure.text(0.09, y, text, fontsize=12.5, weight="bold", color="#176B87", va="top")
            return y - 0.034

        def save_page(pdf, figure):
            pdf.savefig(figure)
            plt.close(figure)

        with PdfPages(output, metadata={"Title": "Topside sounder retrieval: status and validation evidence"}) as pdf:
            # Page 1: report summary before any plots.
            fig = page("Topside sounder retrieval", 1)
            fig.text(0.09, 0.858, "Status, recovery evidence, and truth-independence audit", fontsize=13,
                     color="#176B87", va="top")
            y = heading(fig, "Executive summary", 0.803)
            y = paragraph(fig, "The retrieval code can generate and fit synthetic vertical and oblique topside ionograms. "
                          "This report evaluates one reduced vertical case and asks whether a low ionogram mismatch "
                          "also means that the electron-density profile was recovered. No external observations are used.", y)
            y = heading(fig, "Principal findings", y - 0.029)
            findings = [
                ("Closed-model recovery", f"A true density scale of {true_scale:.2f} was recovered as "
                 f"{scales[closed_index]:.2f} by a six-point grid search (absolute error "
                 f"{abs(scales[closed_index] - true_scale):.2f})."),
                ("Profile-shape mismatch", f"A separately specified altitude perturbation produced a "
                 f"{profile_nrmse:.1%} profile NRMSE over 150–450 km, despite a small ionogram objective."),
                ("No-return failure", f"An empty truth ionogram and an empty candidate scored zero objective "
                 f"while their profiles differed by {no_return_nrmse:.1%} NRMSE."),
            ]
            for label, detail in findings:
                fig.text(0.10, y, label, weight="bold", fontsize=10.5, va="top", color="#173449")
                y = paragraph(fig, detail, y - 0.023, x=0.12, width=85)
                y -= 0.018
            y = heading(fig, "Validation verdict", y - 0.005)
            y = paragraph(fig, "The current demo's truth is generated by the same density parameterization, "
                          "PyIRI background, ray tracer, and image formation used in retrieval. The alternate "
                          "profile test changes only the density formula. Fully independent truth has not "
                          "been established.", y)
            fig.text(0.09, 0.12, "Read pages 2–6 for setup, paired ionograms, interpretation, and the validation gate.",
                     fontsize=9.5, color="#5B7180")
            save_page(pdf, fig)

            # Page 2: methods and an explicit provenance table.
            fig = page("1. Methods and provenance", 2)
            y = paragraph(fig, "One generated 800 km Earth orbit supplies a vertical transmitter/receiver case. "
                          "The background is PyIRI; PHaRLAP traces three frequencies (4, 5, and 6 MHz). "
                          "Only density scale is searched, at 0.80, 0.90, 1.00, 1.10, 1.20, and 1.30. "
                          "This is a grid-search diagnostic, not a run of the full five-parameter optimizer.", 0.855)
            y = heading(fig, "Where each truth comes from", y - 0.032)
            table_ax = fig.add_axes([0.09, 0.48, 0.82, 0.25]); table_ax.axis("off")
            cells = [
                ["Closed model", "Fitted scalar family", "Shared", "Shared", "No"],
                ["Shape mismatch", "Separate altitude formula", "Shared", "Shared", "No"],
                ["No-return control", "Separate altitude formula", "Shared", "Shared", "No"],
            ]
            table = table_ax.table(cellText=cells,
                                   colLabels=["Case", "Density truth", "PyIRI", "PHaRLAP", "External data"],
                                   cellLoc="left", colLoc="left", loc="center",
                                   colWidths=[0.21, 0.34, 0.13, 0.16, 0.13])
            table.auto_set_font_size(False); table.set_fontsize(8.5); table.scale(1, 2.35)
            for (row, _), cell in table.get_celld().items():
                cell.set_edgecolor("#D8E1E6")
                cell.set_facecolor("#DCEAF0" if row == 0 else ("#F5F8FA" if row % 2 else "white"))
                if row == 0:
                    cell.get_text().set_weight("bold")
            fig.text(0.09, 0.468, "Table 1. Truth provenance. “Shared” means the truth generator and candidate model use the same component.",
                     fontsize=8.5, color="#415868", va="top")
            y = heading(fig, "Scoring and accuracy measures", 0.415)
            y = paragraph(fig, "The ionogram objective combines normalized total-power error, O/X split error, "
                          "and Doppler error. The profile NRMSE is root-mean-square density error divided by "
                          "root-mean-square truth density, evaluated at the center grid column from 150 to 450 km. "
                          "The ionogram objective has no direct percent interpretation.", y)
            y = heading(fig, "Independence conclusion", y - 0.02)
            paragraph(fig, "A separate formula is useful for testing shape mismatch, but both sides still "
                      "share Earth geometry, PyIRI, PHaRLAP, and ionogram construction. It is not a "
                      "fully independent ground truth.", y)
            save_page(pdf, fig)

            # Page 3: put the truth and retrieved ionograms next to one another.
            fig = page("2. Closed-model ionograms", 3)
            paragraph(fig, "The synthetic truth uses density scale 1.12. A search over six candidate scales "
                      "selected 1.00. Both ionograms below were re-rendered on a 0.25 MHz display grid "
                      "after selection; only 4, 5, and 6 MHz determined the selected parameter.", 0.855)
            ionogram_pair(fig, closed_truth_display.cases[0].total_image,
                          closed_fit_display.cases[0].total_image)
            y = paragraph(fig, "Figure 1. Synthetic truth (left) and retrieved (right) total-power ionograms "
                          "for the closed-model test. Axes and color limits are identical. Each image is "
                          "normalized to its own peak, so the comparison shows return location and shape, "
                          "not absolute received power.", 0.385, size=9.5, width=98, line_height=0.021,
                          color="#415868")
            y = heading(fig, "What to look for", y - 0.018)
            paragraph(fig, "Compare the frequency and virtual-range locations of the bright returns. "
                      "Their differences explain why the search did not select the generating scale's "
                      "nearest candidate. The numerical objective and parameter error appear on page 4.", y)
            save_page(pdf, fig)

            # Page 4: closed-model objective and parameter recovery.
            fig = page("3. Closed-model recovery", 4)
            paragraph(fig, "First, the truth density is generated by the same scalar family searched by the "
                      "retrieval. The generating scale (1.12) is deliberately absent from the six candidates.", 0.855)
            ax = fig.add_axes([0.15, 0.39, 0.70, 0.36])
            ax.plot(scales, closed_costs, "o-", color="#176B87", linewidth=1.8)
            ax.axvline(true_scale, color="#B6483C", linestyle="--", label=f"Generated {true_scale:.2f}")
            ax.axvline(scales[closed_index], color="#1A774E", linestyle=":", label=f"Selected {scales[closed_index]:.2f}")
            ax.set(xlabel="Candidate density scale", ylabel="Ionogram objective")
            ax.grid(alpha=0.25); ax.legend(fontsize=9)
            y = paragraph(fig, "Figure 2. The six-point search selected scale "
                          f"{scales[closed_index]:.2f} for truth {true_scale:.2f}. "
                          "The objective is not smooth: neighboring candidates can change the available "
                          "ray returns. This plot shows the complete set of evaluated candidates, not a "
                          "continuous fit curve.", 0.344, size=9.5, width=99, line_height=0.021,
                          color="#415868")
            y = heading(fig, "Interpretation", y - 0.016)
            paragraph(fig, f"The absolute scale error is {abs(scales[closed_index] - true_scale):.2f} "
                      f"({100 * abs(scales[closed_index] - true_scale) / true_scale:.1f}% of truth). "
                      "Because the generator and fitter share the same physical and numerical assumptions, "
                      "this is a check of recoverability in one favorable model family, not independent validation.", y)
            save_page(pdf, fig)

            # Page 5: shape mismatch, with density profiles and paired ionograms.
            fig = page("4. Shape mismatch: profile and ionograms", 5)
            paragraph(fig, "Here truth uses an altitude-dependent density factor, 1.10 + 0.10 "
                      "exp[-0.5((h − 320 km)/55 km)²]. Candidates can only multiply the PyIRI background "
                      "by a constant. The search selected scale 1.10.", 0.855)
            ax = fig.add_axes([0.19, 0.60, 0.62, 0.16])
            ax.plot(truth_profile[profile_mask], heights[profile_mask], label="Truth", color="#B6483C", linewidth=1.8)
            ax.plot(fit_profile[profile_mask], heights[profile_mask], label="Retrieved", color="#176B87", linewidth=1.8)
            ax.set(xlabel="Electron density (cm$^{-3}$)", ylabel="Altitude (km)")
            ax.grid(alpha=0.25); ax.legend(fontsize=8)
            paragraph(fig, f"Figure 3. Truth and retrieved density differ by {profile_nrmse:.1%} NRMSE "
                      "over 150–450 km at the center grid column.", 0.548, size=9.5,
                      width=98, line_height=0.021, color="#415868")
            ionogram_pair(fig, shape_truth_display.cases[0].total_image,
                          shape_fit_display.cases[0].total_image, bottom=0.22, height=0.25)
            paragraph(fig, "Figure 4. Synthetic truth (left) and retrieved (right) total-power ionograms, "
                      "re-rendered every 0.25 MHz with common axes and color limits. Each is normalized "
                      "to its own peak. The 4, 5, and 6 MHz retrieval objective is "
                      f"{structural_costs[structural_index]:.5g}, yet the vertical profile still differs "
                      f"by {profile_nrmse:.1%}. The truth formula is separate; the background and ray tracer are shared.",
                      0.173, size=9.5, width=98, line_height=0.021, color="#415868")
            save_page(pdf, fig)

            # Page 6: explain the zero-score failure and state the practical gate.
            fig = page("5. No-return failure and next gate", 6)
            paragraph(fig, "A second altitude-dependent truth gives no returns at 4–6 MHz. The scale-0.80 "
                      "candidate also gives no returns, so both ionograms are empty and their objective is zero.", 0.855)
            ax = fig.add_axes([0.18, 0.48, 0.64, 0.27])
            ax.plot(no_return_profile[profile_mask], heights[profile_mask], label="Generated no-return truth",
                    color="#B6483C", linewidth=1.8)
            ax.plot(zero_fit_profile[profile_mask], heights[profile_mask], label="Scale-0.80 candidate",
                    color="#176B87", linewidth=1.8)
            ax.set(xlabel="Electron density (cm$^{-3}$)", ylabel="Altitude (km)")
            ax.grid(alpha=0.25); ax.legend(fontsize=8)
            y = paragraph(fig, f"Figure 5. Profiles differ by {no_return_nrmse:.1%} NRMSE from 150 to 450 km "
                          "while the ionogram objective is exactly zero. The score cannot distinguish two "
                          "empty observations. Such cases must be flagged and excluded from accuracy claims.",
                          0.429, size=9.5, width=98, line_height=0.021, color="#415868")
            y = heading(fig, "Conclusion", y - 0.014)
            y = paragraph(fig, "The present evidence supports a working synthetic diagnostic and exposes "
                          "weak identifiability. It does not confirm independent truth or measured "
                          "retrieval accuracy.", y)
            y = heading(fig, "Required validation", y - 0.018)
            paragraph(fig, "Use an external density field or measured sounder profiles with matching "
                      "geometry. Keep the truth source outside the retrieval forward model, include "
                      "frequencies with detected returns, and report density bias and NRMSE across "
                      "multiple cases.", y)
            save_page(pdf, fig)

    print(json.dumps(metrics, indent=2))
    print(output)


if __name__ == "__main__":
    main()
