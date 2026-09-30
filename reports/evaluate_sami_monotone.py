"""Evaluate a frozen ionogram fit against the withheld SAMI density field."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "reports")]
from evaluate_general_field_sami import draw_ionogram, metrics, profiles  # noqa: E402
from freeze_sami_monotone_selection import sha256  # noqa: E402

CASE = ROOT / "reports/data/sami3_wave_20170111_0600"
RUN = ROOT / "reports/data/sami3_monotone_topside"
INDICES = [1, 5, 9, 12, 14, 17, 20]


def main() -> None:
    selection = json.loads((RUN / "selection.json").read_text())
    scores = json.loads((RUN / "full_scores.json").read_text())
    plan = json.loads((RUN / "plan.json").read_text())
    if (selection["selection_uses_truth_density"] or
            scores["selection_uses_truth_density"] or
            selection["profile_indices"] != INDICES or
            scores["profile_indices"] != INDICES or
            selection["chosen_candidate"] != scores["winner"] or
            selection["score_file_sha256"] != sha256(RUN / "full_scores.json") or
            selection["plan_file_sha256"] != sha256(RUN / "plan.json")):
        raise ValueError("Selection must be frozen from all seven ionograms")
    name = selection["chosen_candidate"]
    candidate = next(row for row in plan["candidates"] if row["name"] == name)
    selected_path = ROOT / candidate["grid"]
    if selection["selected_grid_sha256"] != sha256(selected_path):
        raise ValueError("Selected grid changed after the ionogram-only choice")
    rows = json.loads((CASE / "manifest.json").read_text())["profiles"]
    locations = np.array([(row["latitude_deg"], row["longitude_deg"])
                          for row in rows])
    alt, truth = profiles(CASE / "truth_grid.nc", locations)
    selected_alt, retrieved = profiles(selected_path, locations)
    if not np.array_equal(alt, selected_alt):
        raise ValueError("Different altitude axes")
    observed = np.array(INDICES) - 1

    insitu = json.loads((RUN / "insitu_1km.json").read_text())
    dense_locations = np.asarray(insitu["lat_lon_deg"])
    _, sampled = profiles(selected_path, dense_locations)
    spacecraft_alt = float(insitu["altitude_km"])
    modeled_insitu = sampled[:, np.where(alt == spacecraft_alt)[0][0]]
    measured_insitu = np.asarray(insitu["density_cm3"])
    insitu_error = 100 * (modeled_insitu / measured_insitu - 1)
    grid_topside_upsteps = int(sum(np.count_nonzero(
        np.diff(row[np.argmax(row):]) > 1e-8 * np.max(row)) for row in retrieved))
    dense_topside_upsteps = int(sum(np.count_nonzero(
        np.diff(row[np.argmax(row):]) > 1e-8 * np.max(row)) for row in sampled))

    result = {
        "selection_uses_truth_density": False,
        "selected_candidate": name,
        "observed_profile_indices": INDICES,
        "model_truth": "SAMI3/HIAMCM electron density remapped onto PyIRI magnetic/collision grid",
        "synthetic_insitu_samples": len(measured_insitu),
        "maximum_insitu_absolute_error_percent": float(np.max(abs(insitu_error))),
        "topside_upward_steps_at_twenty_profile_locations": grid_topside_upsteps,
        "topside_upward_steps_at_2225_insitu_locations": dense_topside_upsteps,
        "observed_positions": metrics(alt, truth[observed], retrieved[observed]),
        "twenty_position_track": metrics(alt, truth, retrieved),
        "ionogram_score": next(row for row in scores["candidates"]
                                if row["name"] == name),
    }
    path = RUN / "evaluation.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    path.chmod(0o600)

    lat = locations[:, 0]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True,
                             sharex=True, sharey=True)
    for ax, field, title in zip(axes[:2], (truth, retrieved),
                                ("SAMI3/HIAMCM truth", "Retrieved density")):
        image = ax.contourf(lat, alt, np.log10(np.maximum(field.T, 1.0)),
                            levels=np.linspace(4.0, 5.7, 18), cmap="viridis",
                            extend="both")
        ax.scatter(lat[observed], np.full(len(observed), 790), s=10,
                   color="white", edgecolors="black", linewidths=.3)
        ax.set(title=title, xlabel="Latitude (degrees)", ylim=(150, 800))
    fig.colorbar(image, ax=axes[:2], label="log10 density (cm$^{-3}$)")
    difference = 100 * (retrieved - truth) / truth.max(axis=1)[:, None]
    image = axes[2].contourf(lat, alt, difference.T, levels=np.linspace(-40, 40, 17),
                              cmap="RdBu_r", extend="both")
    axes[2].set(title="Difference / true peak", xlabel="Latitude (degrees)")
    fig.colorbar(image, ax=axes[2], label="Percent")
    axes[0].set_ylabel("Altitude (km)")
    fig.savefig(RUN / "density_cut.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True,
                             sharex=True, sharey=True)
    for ax, index in zip(axes, (9, 14)):
        row = index - 1
        ax.plot(truth[row] / 1e5, alt, color="black", label="Truth")
        ax.plot(retrieved[row] / 1e5, alt, color="#2166ac", label="Retrieved")
        ax.set(title=f"Profile {index}", xlabel="Density ($10^5$ cm$^{-3}$)",
               ylim=(150, 800))
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Altitude (km)")
    axes[0].legend()
    fig.savefig(RUN / "density_profiles.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.7), constrained_layout=True)
    sounding_lat = lat[observed]
    axes[0].plot(sounding_lat, candidate["fraction_at_quarter_by_profile"],
                 "o-", label="25% of peak-to-spacecraft height")
    axes[0].plot(sounding_lat, candidate["fraction_at_three_fifths_by_profile"],
                 "s-", label="60% of peak-to-spacecraft height")
    axes[0].set(xlabel="Latitude (degrees)",
                ylabel="Fraction of peak-to-spacecraft log-density span")
    axes[0].legend(fontsize=8)
    axes[1].plot(sounding_lat, candidate["peak_height_shift_km_by_profile"],
                 "o-", color="#b2182b")
    axes[1].set(xlabel="Latitude (degrees)",
                ylabel="F2 peak-height shift (km)")
    for ax in axes:
        ax.grid(alpha=.2)
    fig.savefig(RUN / "topside_parameters.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    draw_ionogram(axes[0], CASE / "vertical_truth/ionogram_14.npz",
                   "SAMI3/HIAMCM truth, profile 14")
    draw_ionogram(axes[1], RUN / "forward" / name / "ionograms/ionogram_14.npz",
                   "Retrieved, profile 14")
    axes[0].legend()
    fig.savefig(RUN / "ionogram_14.png", dpi=200)
    plt.close(fig)
    print(json.dumps({"selected": name,
                      "observed": result["observed_positions"],
                      "insitu_max_abs_percent": result["maximum_insitu_absolute_error_percent"],
                      "topside_upsteps_at_2225_locations": dense_topside_upsteps}, indent=2))


if __name__ == "__main__":
    main()
