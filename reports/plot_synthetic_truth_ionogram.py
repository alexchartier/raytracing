"""Plot accepted closed-model truth returns as occupied bins without modeled power.

Run from any directory with ``python3 reports/plot_synthetic_truth_ionogram.py``.
The checked-in NPZ contains the accepted returns from a density-scale-1.12
synthetic truth sweep (2–10 MHz, 100 kHz spacing) with adaptive homing.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "data" / "synthetic_truth_returns_adaptive_2-10MHz_100kHz.npz"
OUTPUT_STEM = HERE / "figures" / "synthetic_truth_ionogram_2-10MHz_100kHz_1km"
FREQUENCIES_MHZ = np.arange(2.0, 10.0001, 0.1)
RANGE_MIN_KM = 150


def main() -> None:
    with np.load(SOURCE, allow_pickle=False) as source:
        records = np.asarray(source["records"], dtype=float)
        method = str(source["method"])
        fan_directions = int(source["fan_launch_directions"])
    assert records.shape[1] == 5
    included = records[records[:, 2] >= RANGE_MIN_KM]
    if included.size == 0:
        raise ValueError("No accepted returns at or above 150 km")
    range_max_km = max(1200, int(np.ceil(included[:, 2].max() / 50.0) * 50 + 50))

    occupied = np.zeros((len(FREQUENCIES_MHZ), range_max_km - RANGE_MIN_KM), dtype=bool)
    for frequency_index, _mode, range_km, _miss_m, _absorption_db in included:
        frequency_index = int(frequency_index)
        if not 0 <= frequency_index < len(FREQUENCIES_MHZ):
            raise ValueError("Accepted return has an invalid frequency index")
        range_index = int(np.floor(range_km)) - RANGE_MIN_KM
        occupied[frequency_index, range_index] = True
    if not np.any(occupied):
        raise ValueError("No accepted returns at or above 150 km")

    fig, ax = plt.subplots(figsize=(12.0, 8.0), facecolor="white")
    fig.subplots_adjust(left=0.10, right=0.96, bottom=0.16, top=0.87)
    ax.set_facecolor("#100b23")
    for frequency_index, range_index in np.argwhere(occupied):
        ax.add_patch(Rectangle(
            (FREQUENCIES_MHZ[frequency_index] - 0.05, RANGE_MIN_KM + range_index),
            0.1, 1.0, facecolor="#ffd166",
            edgecolor="none", antialiased=False,
        ))
    ax.set_xlim(1.95, 10.05)
    ax.set_ylim(range_max_km, RANGE_MIN_KM)
    ax.set_xticks(np.arange(2, 11, 1))
    ax.set_yticks(np.arange(150, range_max_km + 1, 150))
    ax.set_xlabel("Frequency (MHz)", fontsize=13)
    ax.set_ylabel("Virtual range (km)", fontsize=13)
    ax.tick_params(labelsize=11)
    fig.suptitle("Synthetic truth ionogram", x=0.10, y=0.96, ha="left",
                 fontsize=18, fontweight="bold", color="#173449")
    fig.text(0.10, 0.916,
             f"2–10 MHz  ·  100 kHz × 1 km bins  ·  {method} homing, {fan_directions}-direction anchors",
             fontsize=11, color="#435969")
    fig.text(0.10, 0.065,
             f"{len(included)} accepted O/X returns in {np.count_nonzero(occupied)} occupied bins across "
             f"{np.unique(included[:, 0]).size} frequencies. Each occupied bin is shown once; no smoothing.",
             fontsize=10, color="#435969")

    OUTPUT_STEM.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_STEM.with_suffix(".png"), dpi=300)
    fig.savefig(OUTPUT_STEM.with_suffix(".pdf"))
    plt.close(fig)
    print(OUTPUT_STEM.with_suffix(".png"))
    print(OUTPUT_STEM.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
