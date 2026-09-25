"""Plot four vertical ionograms from accepted-return NPZ files.

Example: python3 reports/compare_vertical_fans.py A.npz B.npz C.npz D.npz
Each accepted O/X return is drawn as a 0.1 MHz by 1 km rectangle. No power
estimate, smoothing, or retrieval overlay is used.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch, Rectangle
from scipy.optimize import linear_sum_assignment


def read_sweep(path: Path) -> dict:
    with np.load(path, allow_pickle=False) as source:
        records = np.asarray(source["records"], dtype=float)
        frequencies = np.asarray(source["frequencies_mhz"], dtype=float)
        return {
            "records": records[records[:, 2] >= 150.0],
            "frequencies": frequencies,
            "fan": int(source["fan_launch_directions"]),
            "anchor_fan": int(source["anchor_fan_launch_directions"]),
            "runtime": float(source["runtime_seconds"]) if "runtime_seconds" in source else None,
            "tolerance": float(source["homing_tolerance_m"]),
            "layout": str(source["vertical_fan_layout"]),
            "outer_fraction": float(source["vertical_outer_ray_fraction"]),
            "guard_seed_limit": int(source["vertical_guard_seed_limit"]),
        }


def match_count(reference: np.ndarray, candidate: np.ndarray) -> int:
    """Count one-to-one O/X returns within 1 km at the same frequency."""
    matched = 0
    for frequency_index in range(81):
        for mode in (-1, 1):
            first = reference[(reference[:, 0] == frequency_index) & (reference[:, 1] == mode), 2]
            second = candidate[(candidate[:, 0] == frequency_index) & (candidate[:, 1] == mode), 2]
            if first.size and second.size:
                distances = np.abs(first[:, None] - second[None, :])
                rows, columns = linear_sum_assignment(np.where(distances <= 1.0, distances, 1e6))
                matched += int(np.count_nonzero(distances[rows, columns] <= 1.0))
    return matched


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sweeps", type=Path, nargs=4, help="A, B, C, and D NPZ files")
    parser.add_argument("--output-stem", type=Path,
                        default=Path(__file__).resolve().parent / "figures" / "vertical_fan_abcd")
    args = parser.parse_args()
    sweeps = [read_sweep(path) for path in args.sweeps]
    for sweep in sweeps:
        np.testing.assert_allclose(sweep["frequencies"], np.arange(2.0, 10.0001, 0.1))
        if sweep["tolerance"] != 1000.0:
            raise ValueError("Comparison expects the unchanged 1,000 m homing gate")
    for index, sweep in enumerate(sweeps):
        expected_layout = "az_el" if index == 0 else "equal_area_guarded"
        expected_fraction = (1.0, 1.0, 0.5, 0.5)[index]
        expected_seed_limit = (-1, 0, 0, 4)[index]
        if (sweep["layout"] != expected_layout or
                sweep["outer_fraction"] != expected_fraction or
                sweep["guard_seed_limit"] != expected_seed_limit):
            raise ValueError(f"Sweep {index + 1} has the wrong fan settings")
    reference = sweeps[0]["records"]
    maximum = max(float(sweep["records"][:, 2].max()) for sweep in sweeps)
    range_max = int(np.ceil(maximum / 50.0) * 50.0 + 50.0)

    fig, axes = plt.subplots(2, 2, figsize=(19, 14), sharex=True, sharey=True)
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.09, top=0.865, wspace=0.14, hspace=0.25)
    colors = {1: "#ffd166", -1: "#ef7d67"}
    for index, (ax, sweep) in enumerate(zip(axes.flat, sweeps)):
        ax.set_facecolor("#100b23")
        records = sweep["records"]
        for frequency_index, mode, range_km, _miss_m, _absorption_db in records:
            # Draw every accepted record, including coincident O/X returns.
            ax.add_patch(Rectangle(
                (float(sweep["frequencies"][int(frequency_index)]) - 0.05, np.floor(range_km)),
                0.1, 1.0, facecolor=colors[int(mode)], edgecolor="none", alpha=0.8,
                antialiased=False,
            ))
        ax.set_xlim(1.95, 10.05)
        ax.set_ylim(range_max, 150.0)
        ax.set_xticks(np.arange(2, 11, 1))
        ax.set_yticks(np.arange(150, range_max + 1, 150))
        ax.tick_params(labelsize=12)
        letter = "ABCD"[index]
        runtime = "runtime unavailable" if sweep["runtime"] is None else f'{sweep["runtime"]:.1f} s'
        matched = match_count(reference, records)
        comparison = "reference" if index == 0 else f"{matched}/{len(reference)} reference returns matched"
        name = ("Original az/el", "Guarded full outer fan, 0 extra seeds",
                "Guarded half outer fan, 0 extra seeds",
                "Guarded half outer fan, 4 extra seeds")[index]
        ax.set_title(
            f"{letter}  {name}  ·  {len(records)} returns  ·  {runtime}\n"
            f'{sweep["anchor_fan"]} anchor / {sweep["fan"]} full directions  ·  {comparison}',
            loc="left", fontsize=14, color="#173449", pad=11,
        )
        if index // 2 == 1:
            ax.set_xlabel("Frequency (MHz)", fontsize=13)
        if index % 2 == 0:
            ax.set_ylabel("Virtual range (km)", fontsize=13)

    fig.suptitle("Vertical ionogram fan comparison", x=0.075, y=0.965,
                 ha="left", fontsize=24, fontweight="bold", color="#173449")
    fig.text(0.075, 0.921,
             "Synthetic ionosphere · 2–10 MHz in 100 kHz steps · 1 km range bins · all accepted O/X returns",
             fontsize=14, color="#435969")
    fig.legend(handles=[Patch(facecolor=colors[1], label="O mode"),
                        Patch(facecolor=colors[-1], label="X mode")],
               loc="upper right", bbox_to_anchor=(0.985, 0.945), ncol=2, frameon=False, fontsize=13)
    fig.text(0.075, 0.035,
             "Homing gate: 1,000 m throughout. Match: same frequency and mode, virtual ranges within 1 km; "
             "coincident returns may blend. Timings include setup and both modes.",
             fontsize=11, color="#435969")
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_stem.with_suffix(".png"), dpi=600)
    fig.savefig(args.output_stem.with_suffix(".pdf"))
    plt.close(fig)
    print(args.output_stem.with_suffix(".png"))
    print(args.output_stem.with_suffix(".pdf"))
    for letter, sweep in zip("ABCD", sweeps):
        print(f'{letter}: {len(sweep["records"])} returns, '
              f'{match_count(reference, sweep["records"])} matched, '
              f'{sweep["runtime"]} seconds, '
              f'{sweep["anchor_fan"]}/{sweep["fan"]} anchor/full directions')


if __name__ == "__main__":
    main()
