"""Plot truth and retrieved D ionograms side by side with O/X modes overlaid.

Every accepted return is drawn in its 0.1 MHz by 1 km display cell. Repeated
returns in a cell darken that cell rather than being discarded.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import numpy as np

from fit_d_ionogram import Ionogram


COLORS = {1: "#2455a4", -1: "#c43b31"}


def plot_pair(truth_path: Path, retrieved_path: Path, output: Path, *,
              truth_label: str = "Truth", retrieved_label: str = "Retrieved") -> None:
    truth = Ionogram.read(truth_path)
    retrieved = Ionogram.read(retrieved_path)
    if not np.array_equal(truth.frequencies, retrieved.frequencies):
        raise ValueError("Frequency axes differ")

    maximum = max(float(np.max(ionogram.records[:, 2]))
                  for ionogram in (truth, retrieved))
    range_max_km = max(1600, int(np.ceil((maximum + 1) / 100) * 100))
    fig, axes = plt.subplots(1, 2, figsize=(15, 11), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, ionogram, title in zip(
            axes, (truth, retrieved), (truth_label, retrieved_label)):
        # Keep one patch per record so duplicate accepted returns remain in the
        # image. Alpha makes cells with multiple returns visibly darker.
        for record in ionogram.records:
            index, mode, group_range_km = int(record[0]), int(record[1]), float(record[2])
            if mode not in COLORS:
                raise ValueError(f"Unexpected propagation mode: {mode}")
            frequency = float(ionogram.frequencies[index])
            ax.add_patch(Rectangle(
                (frequency - 0.05, np.floor(group_range_km)), 0.1, 1.0,
                facecolor=COLORS[mode], edgecolor="none", alpha=0.82,
            ))
        ax.set(title=f"{title} ({len(ionogram.records)} accepted returns)",
               xlabel="Frequency (MHz)", xlim=(1.95, 10.05), ylim=(150, range_max_km))
        ax.grid(alpha=0.15, linewidth=0.5)
    axes[0].set_ylabel("Group range (km)")
    axes[0].legend(handles=[Patch(facecolor=COLORS[1], label="O mode (blue)"),
                            Patch(facecolor=COLORS[-1], label="X mode (red)")],
                   loc="upper right", frameon=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("truth", type=Path)
    parser.add_argument("retrieved", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--truth-label", default="Truth")
    parser.add_argument("--retrieved-label", default="Retrieved")
    args = parser.parse_args()
    plot_pair(args.truth, args.retrieved, args.output,
              truth_label=args.truth_label, retrieved_label=args.retrieved_label)


if __name__ == "__main__":
    main()
