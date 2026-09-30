"""Evaluate the frozen two-satellite fit at link midpoints and make figures."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "reports/data/sami3_wave_20170111_0600"
RUN = ROOT / "reports/data/sami3_monotone_oblique"
os.environ["SOUNDER_CASE_ROOT"] = str(CASE)
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from evaluate_general_field_sami import metrics, profiles  # noqa: E402
from freeze_sami_monotone_selection import sha256  # noqa: E402
from oblique_wave_pass import endpoints  # noqa: E402
from python_raytrace.grid import load_ionosphere_grid_netcdf  # noqa: E402

INDICES = [1, 5, 9, 12, 14, 17, 20]
OLD_GRID = CASE / "oblique_wave_600km/gain1_sigma60_height0.15_spline/grid.nc"
TRANSFER_GRID = ROOT / "reports/data/sami3_monotone_topside/varying_g04/grid.nc"


def draw_ionogram(ax, path: Path, title: str) -> None:
    with np.load(path, allow_pickle=False) as saved:
        frequency = saved["frequencies_mhz"]
        returns = saved["records"]
    low, high = 1000, 1900
    extent = (float(frequency[0] - .05), float(frequency[-1] + .05),
              low - .5, high + .5)
    for mode, color, label in ((1, "#2166ac", "O"), (-1, "#b2182b", "X")):
        selected = returns[returns[:, 1] == mode]
        pixels = np.zeros((high - low + 1, len(frequency), 4), dtype=np.float32)
        ranges = np.rint(selected[:, 2]).astype(int)
        indices = selected[:, 0].astype(int)
        valid = (ranges >= low) & (ranges <= high)
        pixels[ranges[valid] - low, indices[valid]] = to_rgba(color, alpha=.8)
        ax.imshow(pixels, origin="lower", extent=extent, interpolation="nearest",
                  aspect="auto", label=label)
    ax.set(xlim=(2, 10), ylim=(low, high), xlabel="Frequency (MHz)",
           ylabel="Group range (km)", title=title)
    ax.grid(alpha=.2)


def upward_steps(field: np.ndarray) -> int:
    count = 0
    for row in field:
        peak = int(np.argmax(row))
        count += int(np.count_nonzero(np.diff(row[peak:]) > 1e-8 * row[peak]))
    return count


def main() -> None:
    os.umask(0o077)
    selection = json.loads((RUN / "selection.json").read_text())
    scores = json.loads((RUN / "full_scores.json").read_text())
    feasible = [row for row in scores["candidates"]
                if selection["insitu_max_absolute_error_percent_by_candidate"][row["name"]]
                <= selection["insitu_hard_tolerance_percent"]]
    if (selection["selection_uses_truth_density"] or
            selection["profile_indices"] != INDICES or
            not feasible or
            selection["chosen_candidate"] != min(
                feasible, key=lambda row: row["mean_total_equivalent_km"])["name"] or
            selection["score_file_sha256"] != sha256(RUN / "full_scores.json") or
            selection["plan_file_sha256"] != sha256(RUN / "plan.json")):
        raise ValueError("Oblique ionogram selection must be frozen first")
    chosen_grid = ROOT / selection["grid"]
    if selection["grid_sha256"] != sha256(chosen_grid):
        raise ValueError("Selected grid changed after selection")

    all_links = np.array([[(endpoints(i)[0].lat_deg + endpoints(i)[1].lat_deg) / 2,
                           endpoints(i)[0].lon_deg] for i in range(1, 21)])
    selected_links = all_links[np.array(INDICES) - 1]
    altitude, truth = profiles(CASE / "truth_grid.nc", all_links)
    _, old = profiles(OLD_GRID, all_links)
    _, transfer = profiles(TRANSFER_GRID, all_links)
    _, retrieved = profiles(chosen_grid, all_links)
    observed = np.array(INDICES) - 1

    samples = json.loads((RUN / "insitu_1km.json").read_text())
    _, track = profiles(chosen_grid, np.asarray(samples["lat_lon_deg"]))
    at_spacecraft = np.array([np.interp(samples["altitude_km"], altitude, row)
                              for row in track])
    track_error = 100 * (at_spacecraft / np.asarray(samples["density_cm3"]) - 1)
    grid = load_ionosphere_grid_netcdf(chosen_grid)
    full_upsteps = upward_steps(grid.iono_en_grid.reshape(-1, len(altitude)))
    result = {
        "geometry": "two spacecraft at 800 km, 600 km along-track separation",
        "ionogram_content": "O/X reflected returns only; quasi-direct path excluded by generator",
        "frequency_grid": "2-10 MHz at 100 kHz",
        "selection_uses_truth_density": False,
        "selection_uses_synthetic_insitu_density_samples": True,
        "selection_uses_full_3d_sami_density": False,
        "truth_density_provenance": "SAMI3/HIAMCM electron density, independent of the PyIRI prior; shared PyIRI magnetic and collision fields",
        "development_is_blind": False,
        "chosen_candidate": selection["chosen_candidate"],
        "selected_grid": selection["grid"],
        "ionogram_scores": {row["name"]: {
            "mean_total_equivalent_km": row["mean_total_equivalent_km"],
            "mean_doppler_1hz_score": row["mean_doppler_1hz_score"]}
            for row in scores["candidates"]},
        "seven_observed_link_midpoints": {
            "old_oblique_fit": metrics(altitude, truth[observed], old[observed]),
            "vertical_transfer": metrics(altitude, truth[observed], transfer[observed]),
            "chosen": metrics(altitude, truth[observed], retrieved[observed]),
        },
        "twenty_link_midpoints": {
            "old_oblique_fit": metrics(altitude, truth, old),
            "vertical_transfer": metrics(altitude, truth, transfer),
            "chosen": metrics(altitude, truth, retrieved),
        },
        "in_situ_sample_count": len(at_spacecraft),
        "in_situ_max_absolute_error_percent": float(np.max(abs(track_error))),
        "topside_upward_steps_in_forward_grid": full_upsteps,
        "topside_upward_steps_at_link_midpoints": upward_steps(retrieved),
    }
    output = RUN / "evaluation.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    output.chmod(0o600)

    latitudes = all_links[:, 0]
    figure, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True,
                                sharex=True, sharey=True)
    for ax, field, title in zip(axes[:2], (truth, retrieved),
                                ("SAMI3/HIAMCM truth", "Selected retrieval")):
        image = ax.contourf(latitudes, altitude, np.log10(np.maximum(field.T, 1.0)),
                            levels=np.linspace(4.0, 5.7, 18), cmap="viridis",
                            extend="both")
        ax.scatter(selected_links[:, 0], np.full(len(INDICES), 790), s=12,
                   color="white", edgecolors="black", linewidths=.3)
        ax.set(title=title, xlabel="Link-midpoint latitude (degrees)", ylim=(150, 800))
    figure.colorbar(image, ax=axes[:2], label="log10 density (cm$^{-3}$)")
    difference = 100 * (retrieved - truth) / truth.max(axis=1)[:, None]
    diff = axes[2].contourf(latitudes, altitude, difference.T,
                            levels=np.linspace(-40, 40, 17), cmap="RdBu_r",
                            extend="both")
    axes[2].set(title="Difference / true peak", xlabel="Link-midpoint latitude (degrees)")
    axes[0].set_ylabel("Altitude (km)")
    figure.colorbar(diff, ax=axes[2], label="Percent")
    figure.savefig(RUN / "density_cut.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    draw_ionogram(axes[0], CASE / "oblique_wave_600km/truth/ionogram_14.npz",
                   "Truth, link 14")
    name = selection["chosen_candidate"]
    modeled = (CASE / "oblique_wave_600km/baseline" if name == "old_baseline" else
               CASE / "oblique_wave_600km/gain1_sigma60_height0.15_spline/ionograms"
               if name == "old_fit" else RUN / "forward" / name / "ionograms")
    draw_ionogram(axes[1], modeled / "ionogram_14.npz", "Retrieved, link 14")
    axes[0].legend(handles=[Patch(color="#2166ac", label="O"),
                            Patch(color="#b2182b", label="X")])
    figure.savefig(RUN / "ionogram_14.png", dpi=300)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(9, 5), constrained_layout=True,
                                sharex=True, sharey=True)
    for ax, index in zip(axes, (9, 14)):
        row = index - 1
        ax.plot(truth[row] / 1e5, altitude, color="black", label="SAMI truth")
        ax.plot(retrieved[row] / 1e5, altitude, color="#2166ac", label="Retrieved")
        ax.set(title=f"Link {index} midpoint", xlabel="Density ($10^5$ cm$^{-3}$)",
               ylim=(150, 800))
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Altitude (km)")
    axes[0].legend()
    figure.savefig(RUN / "density_profiles.png", dpi=180)
    plt.close(figure)
    print(json.dumps({"selection": selection["chosen_candidate"],
                      "selected_scores": result["ionogram_scores"][name],
                      "seven_midpoint_metrics": result["seven_observed_link_midpoints"]["chosen"],
                      "in_situ_max_absolute_error_percent": result["in_situ_max_absolute_error_percent"],
                      "topside_upward_steps_in_forward_grid": full_upsteps}, indent=2))


if __name__ == "__main__":
    main()
