"""Reproduce a small, clearly labeled topside retrieval diagnostic PDF.

Run from the repository root with ``python3 reports/build_topside_validation_report.py``.
Requires the installed PyIRI/PHaRLAP runtime. No AMPERE data file is required.
"""

from __future__ import annotations

import datetime as dt
import hashlib
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
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

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
from python_raytrace import multisat_topside_inverse_demo as inverse_demo  # noqa: E402
from python_raytrace.tracer import PointToPointRayTracer  # noqa: E402


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
        # Plot every accepted O/X return from 2 to 10 MHz. Store the raw return
        # table locally so later plot edits do not need another raytrace sweep.
        display_frequencies = np.arange(2.0, 10.0001, 0.1)
        display_range_edges = np.arange(0.0, config.range_max_km + 1.0, 1.0)
        display_problem = replace(
            problem,
            frequencies_mhz=display_frequencies,
            range_edges_km=display_range_edges,
            range_centers_km=0.5 * (display_range_edges[:-1] + display_range_edges[1:]),
        )

        def render_all_returns(label: str, grid):
            cache_hash = hashlib.sha256()
            cache_hash.update(Path(inverse_demo.__file__).read_bytes())
            for values in (
                grid.iono_en_grid, grid.collision_freq, grid.Bx, grid.By, grid.Bz,
                case.fan_elevations_deg, case.fan_bearings_deg, display_frequencies,
            ):
                cache_hash.update(np.asarray(values, dtype=np.float64).tobytes())
            cache_hash.update(repr((case.tx_points[1], case.rx_points[1],
                                    config.homing_tolerance_m, config.homed_max_returns_per_frequency,
                                    config.seed_max_candidates_per_frequency)).encode())
            cache_dir = ROOT / ".cache" / "topside_report_returns"
            cache_dir.mkdir(parents=True, exist_ok=True)
            cache_path = cache_dir / f"{label.replace(' ', '_')}_{cache_hash.hexdigest()[:16]}.npz"
            if cache_path.exists():
                with np.load(cache_path, allow_pickle=False) as data:
                    records = np.asarray(data["records"], dtype=float)
                    count_array = np.asarray(data["count_array"], dtype=int)
                print(f"Loaded accepted returns: {label}", flush=True)
            else:
                record_list: list[tuple[float, float, float, float, float]] = []
                count_array = np.zeros((display_frequencies.size, 2), dtype=int)
                tracer = PointToPointRayTracer()
                print(f"Tracing all accepted returns, 2–10 MHz at 100 kHz: {label}", flush=True)
                for frequency_index, frequency_mhz in enumerate(display_frequencies):
                    for mode_index, mode in enumerate((1, -1)):
                        returns = inverse_demo._home_frequency_returns(
                            tracer,
                            tx=case.tx_points[1],
                            rx=case.rx_points[1],
                            grid=grid,
                            fan_elevations_deg=case.fan_elevations_deg,
                            fan_bearings_deg=case.fan_bearings_deg,
                            frequency_mhz=float(frequency_mhz),
                            ox_mode=mode,
                            config=config,
                        )
                        count_array[frequency_index, mode_index] = len(returns)
                        record_list.extend(
                            (float(frequency_index), float(mode), ray.group_range_km,
                             ray.miss_m, ray.absorption_db)
                            for ray in returns
                        )
                    if frequency_index % 20 == 0:
                        print(f"  {label}: {frequency_index + 1}/{display_frequencies.size} frequencies", flush=True)
                records = np.asarray(record_list, dtype=float).reshape(-1, 5)
                np.savez_compressed(cache_path, records=records, count_array=count_array)

            counts = {
                (round(float(frequency_mhz), 5), mode): int(count_array[frequency_index, mode_index])
                for frequency_index, frequency_mhz in enumerate(display_frequencies)
                for mode_index, mode in enumerate((1, -1))
            }
            image = np.zeros((display_frequencies.size, display_problem.range_centers_km.size), dtype=float)
            plotted_returns = 0
            for frequency_index, _mode, group_range_km, miss_m, absorption_db in records:
                if not display_problem.range_edges_km[0] <= group_range_km < display_problem.range_edges_km[-1]:
                    raise ValueError(f"Accepted {label} return at {group_range_km:.3f} km falls outside the display grid")
                range_index = int(np.searchsorted(display_problem.range_edges_km, group_range_km, side="right") - 1)
                weight = math.exp(-0.5 * (miss_m / 35_000.0) ** 2) * 10.0 ** (-max(absorption_db, 0.0) / 10.0)
                image[int(frequency_index), range_index] += weight
                plotted_returns += 1
            occupied_bins = int(np.count_nonzero(image))
            peak = float(np.max(image))
            if peak > 0.0:
                image /= peak
            in_retrieval_range = (
                (records[:, 2] >= config.range_min_km)
                & (records[:, 2] < config.range_max_km)
            )
            return image, counts, {
                "all_accepted": int(records.shape[0]),
                "plotted_returns": plotted_returns,
                "occupied_bins": occupied_bins,
                "under_10_km": int(np.count_nonzero(records[:, 2] < 10.0)),
                "below_retrieval_range_floor": int(np.count_nonzero(records[:, 2] < config.range_min_km)),
                "within_retrieval_range": int(np.count_nonzero(in_retrieval_range)),
            }

        closed_truth_grid = _apply_fit_params_to_grid(problem, background, replace(zero_wave, density_scale=true_scale))
        closed_fit_grid = _apply_fit_params_to_grid(
            problem, background, replace(zero_wave, density_scale=float(scales[closed_index]))
        )
        shape_fit_grid = _apply_fit_params_to_grid(
            problem, background, replace(zero_wave, density_scale=float(scales[structural_index]))
        )
        closed_truth_field, closed_truth_counts, closed_truth_return_stats = render_all_returns(
            "closed-model truth", closed_truth_grid
        )
        closed_fit_field, closed_fit_counts, closed_fit_return_stats = render_all_returns(
            "closed-model retrieved", closed_fit_grid
        )
        shape_truth_field, shape_truth_counts, shape_truth_return_stats = render_all_returns(
            "profile-shape truth", independent_grid
        )
        shape_fit_field, shape_fit_counts, shape_fit_return_stats = render_all_returns(
            "profile-shape retrieved", shape_fit_grid
        )
        display_fields = (closed_truth_field, closed_fit_field, shape_truth_field, shape_fit_field)

        def count_summary(counts: dict[tuple[float, int], int]) -> dict[str, int]:
            return {
                "max_returns_per_frequency_and_mode": max(counts.values(), default=0),
                "frequency_mode_cells_with_multiple_returns": sum(value > 1 for value in counts.values()),
                "frequency_mode_cells_with_any_return": sum(value > 0 for value in counts.values()),
            }

        metrics["display"] = {
            "frequency_start_mhz": 2.0,
            "frequency_stop_mhz": 10.0,
            "frequency_step_khz": 100.0,
            "range_bin_km": 1.0,
            "range_start_km": 0.0,
            "retrieval_objective_range_floor_km": config.range_min_km,
            "homing_tolerance_m": config.homing_tolerance_m,
            "max_accepted_returns_per_frequency_and_mode": config.homed_max_returns_per_frequency,
            "display_selection": "all accepted O/X returns in 1 km bins from 0 to 2600 km; coincident returns summed; no smoothing",
            "accepted_return_counts": {
                "closed_truth": closed_truth_return_stats,
                "closed_fit": closed_fit_return_stats,
                "shape_truth": shape_truth_return_stats,
                "shape_fit": shape_fit_return_stats,
            },
            "closed_truth_return_counts": count_summary(closed_truth_counts),
            "closed_fit_return_counts": count_summary(closed_fit_counts),
            "shape_truth_return_counts": count_summary(shape_truth_counts),
            "shape_fit_return_counts": count_summary(shape_fit_counts),
        }
        metrics_path.write_text(json.dumps(metrics, indent=2) + "\n")
        visible_range_indices = np.flatnonzero(np.max(np.stack(display_fields), axis=(0, 1)) > 0.0)
        if visible_range_indices.size:
            display_range_min = max(float(display_problem.range_edges_km[0]),
                                    float(display_problem.range_edges_km[visible_range_indices[0]]) - 100.0)
            display_range_max = min(float(display_problem.range_edges_km[-1]),
                                    float(display_problem.range_edges_km[visible_range_indices[-1] + 1]) + 100.0)
        else:
            display_range_min = float(display_problem.range_edges_km[0])
            display_range_max = float(display_problem.range_edges_km[-1])

        def ionogram_pair(figure, truth_field, fit_field, *, bottom=0.43, height=0.30):
            axes = (figure.add_axes([0.11, bottom, 0.32, height]),
                    figure.add_axes([0.55, bottom, 0.32, height]))
            colormap = plt.get_cmap("magma")
            for ax, field, title in zip(axes, (truth_field, fit_field),
                                        ("Synthetic truth ionogram", "Retrieved ionogram")):
                ax.set_facecolor(colormap(0.0))
                # Keep each occupied 0.1 MHz × 1 km bin as a vector rectangle.
                # A rasterized full-height image drops subpixel near-range bins.
                for frequency_index, range_index in np.argwhere(field > 0.0):
                    ax.add_patch(Rectangle(
                        (display_frequencies[frequency_index] - 0.05,
                         display_problem.range_edges_km[range_index]),
                        0.1, 1.0, facecolor=colormap(field[frequency_index, range_index]),
                        edgecolor="none", antialiased=False,
                    ))
                ax.set_ylim(display_range_max, display_range_min)
                ax.set_xlim(1.95, 10.05)
                ax.axhline(config.range_min_km, color="cyan", linestyle="--", linewidth=0.8, alpha=0.9)
                ax.text(9.95, config.range_min_km + 12.0, "150 km retrieval floor",
                        ha="right", va="top", fontsize=6.5, color="cyan",
                        bbox={"facecolor": "#21122e", "edgecolor": "none", "alpha": 0.8, "pad": 1.5})
                ax.set_title(title, fontsize=10, weight="bold")
                ax.set_xlabel("Frequency (MHz)")
                ax.set_ylabel("Virtual range (km)")
                near_range = ax.inset_axes([0.52, 0.32, 0.45, 0.24])
                near_range.imshow(
                    field[:, :10].T, origin="lower", aspect="auto", cmap=colormap,
                    vmin=0.0, vmax=1.0, extent=(1.95, 10.05, 0.0, 10.0),
                    interpolation="nearest",
                )
                near_range.set_xlim(1.95, 10.05)
                near_range.set_ylim(10.0, 0.0)
                near_range.set_xticks([2, 6, 10])
                near_range.set_yticks([0, 5, 10])
                near_range.tick_params(labelsize=6, colors="cyan", length=2)
                near_range.set_title("0–10 km detail", fontsize=6.5, color="cyan", loc="left", pad=1)
                for spine in near_range.spines.values():
                    spine.set_color("cyan")
                    spine.set_linewidth(0.6)
            colorbar_ax = figure.add_axes([0.89, bottom, 0.015, height])
            figure.colorbar(ScalarMappable(norm=Normalize(0.0, 1.0), cmap=colormap),
                            cax=colorbar_ax, label="Normalized power").ax.tick_params(labelsize=8)

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
                          "The background is PyIRI; the scale search uses 4, 5, and 6 MHz. Only density "
                          "scale is searched, at 0.80 through 1.30 in 0.10 steps. The displayed ionograms "
                          "are re-rendered from 2 to 10 MHz every 100 kHz in 1 km range bins. Homing "
                          "tolerance is 1,000 m, with up to 10 accepted returns per frequency and mode. "
                          "Every accepted O/X return is plotted from 0 km, including short-range returns "
                          "below the 150 km retrieval floor. Returns sharing a bin are summed; no smoothing "
                          "is applied. The retrieval objective uses 150–2,600 km. "
                          "This is a grid-search diagnostic, not the full five-parameter optimizer.", 0.855)
            y = heading(fig, "Where each truth comes from", y - 0.032)
            table_ax = fig.add_axes([0.09, 0.385, 0.82, 0.25]); table_ax.axis("off")
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
            fig.text(0.09, 0.373, "Table 1. Truth provenance. “Shared” means the truth generator and candidate model use the same component.",
                     fontsize=8.5, color="#415868", va="top")
            y = heading(fig, "Scoring and accuracy measures", 0.32)
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
                      "selected 1.00. Both ionograms below were re-rendered from 2 to 10 MHz every "
                      "100 kHz in 1 km range bins; only 4, 5, and 6 MHz selected the parameter.", 0.855)
            ionogram_pair(fig, closed_truth_field, closed_fit_field)
            y = paragraph(fig, "Figure 1. Synthetic truth (left) and retrieved (right) ionograms show every "
                          "accepted O/X return, including those below the dashed 150 km retrieval floor. "
                          "Insets enlarge 0–10 km. Returns in the same 1 km bin are summed, with no smoothing. "
                          "Axes and color limits "
                          "match; each image is normalized to its own "
                          "peak, so absolute received power cannot be compared.", 0.385,
                          size=9.5, width=98, line_height=0.021,
                          color="#415868")
            y = heading(fig, "What to look for", y - 0.018)
            paragraph(fig, "Compare the frequency and virtual-range locations of the bright returns. "
                      f"Of {closed_truth_return_stats['all_accepted']} accepted truth returns, "
                      f"{closed_truth_return_stats['under_10_km']} lie below 10 km and are outside the score. "
                      f"There are {count_summary(closed_truth_counts)['frequency_mode_cells_with_multiple_returns']} "
                      "frequency/mode cells with multiple accepted returns. The 1,000 m homing gate "
                      "checks receiver miss distance; it does not suppress multipath. The parameter "
                      "error and search objective appear on page 4.", y)
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
            ionogram_pair(fig, shape_truth_field, shape_fit_field, bottom=0.22, height=0.25)
            paragraph(fig, "Figure 4. Truth (left) and retrieved (right) ionograms from 2 to 10 MHz at "
                      "100 kHz × 1 km show all accepted O/X returns, including those below the dashed "
                      "150 km retrieval floor; insets enlarge 0–10 km. Coincident returns are summed without smoothing. Each panel "
                      "is normalized to its own peak. The "
                      "4, 5, and 6 MHz search objective is "
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
